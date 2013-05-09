package org.renci.sequencing.util;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;

public class SortedBaseHolder
{
	private LogisticRegression tRegression;
	private StringBuffer tRecalScoresBuffer;
	private List<Integer> tRecalScoresList;
	private Set<String> tSetOfBadBases;
	private boolean bApplyPhredScale = false;

	private int PHRED_SCALE_VALUE = 33;

	private int RECORDS_TO_READ = 1000;
	private BufferedWriter tDebugWriter;
	private LinkedList<SortedBase> localregion;
	private int iCurrentPos = 0;
	private static final int MASK = 255;

	public SortedBaseHolder(List<Float> tCoeffListIn, List<Integer> tFlaggedListIn, Set<String> tSetOfBadBasesIn, boolean bScalePhredBy33In)
			throws IOException
			{
		this.localregion = new LinkedList();
		this.tRegression = LogisticRegression.getInstance(tCoeffListIn, tFlaggedListIn, this.bApplyPhredScale);
		this.tRecalScoresBuffer = new StringBuffer();
		this.tRecalScoresList = new ArrayList();
		this.tSetOfBadBases = tSetOfBadBasesIn;
		this.bApplyPhredScale = bScalePhredBy33In;
			}

	public SortedBaseHolder(List<Float> tCoeffListIn, List<Integer> tFlaggedListIn, Set<String> tSetOfBadBasesIn, boolean bScalePhredBy33In, BufferedWriter tDebugWriterIn) throws IOException
	{
		this.localregion = new LinkedList();
		this.tRegression = LogisticRegression.getInstance(tCoeffListIn, tFlaggedListIn, this.bApplyPhredScale, tDebugWriterIn);
		this.tRecalScoresBuffer = new StringBuffer();
		this.tRecalScoresList = new ArrayList();
		this.tSetOfBadBases = tSetOfBadBasesIn;
		this.bApplyPhredScale = bScalePhredBy33In;
		this.tDebugWriter = tDebugWriterIn;
	}

	public void flushBefore(String chromosome, int location)
			throws IOException
			{
		boolean done = false;
		while ((this.localregion.size() > 0) && (!done))
		{
			SortedBase w;
			if (!((SortedBase)this.localregion.getFirst()).getChromosome().equals(chromosome))
			{
				w = (SortedBase)this.localregion.pop();
			}
			else
			{
				//  SortedBase w;
				if (((SortedBase)this.localregion.getFirst()).getRefLocation() < location)
				{
					w = (SortedBase)this.localregion.pop();
				}
				else
				{
					done = true;
				}
			}
		}
			}

	public void doRegression(String sBaseString, SortedBase tSortedBaseCurrent, int iTrackIt, char tTestBaseChar, int[] iArrayOfQualityScores, double dAverageQualityScore, boolean bNegativeStrand, int iReadLength, Iterator<SortedBase> tSortedBaseIterator, double dNewPhredScore, double dOldPhredScore, RegressionScoreHolder tHolder, List<RegressionScoreHolder> tHolderList)
			throws Exception
			{
		boolean bIsBadBase = false;

		tSortedBaseCurrent.add(iTrackIt, tTestBaseChar, iArrayOfQualityScores[iTrackIt], dAverageQualityScore, bNegativeStrand, iReadLength);
		if (tSortedBaseIterator.hasNext())
		{
			tSortedBaseCurrent = (SortedBase)tSortedBaseIterator.next();

			if (isGoodBase(tTestBaseChar))
			{
				int iNewReadPos = iTrackIt + 1;

				this.tRegression.calcRowVector(iArrayOfQualityScores[iTrackIt], dAverageQualityScore, tTestBaseChar, iNewReadPos, 0, bNegativeStrand, iReadLength);

				dNewPhredScore = this.tRegression.runRegression();

				tHolder = RegressionScoreHolder.makeHolder(Character.toString(tTestBaseChar), new Float((float)dNewPhredScore), sBaseString, bIsBadBase, null);

				tHolderList.add(tHolder);
			}
			else
			{
				bIsBadBase = true;

				dNewPhredScore = iArrayOfQualityScores[iTrackIt];

				tHolder = RegressionScoreHolder.makeHolder(Character.toString(sBaseString.charAt(iTrackIt)), new Float((float)dOldPhredScore), sBaseString, bIsBadBase, new Integer(iArrayOfQualityScores[iTrackIt]));

				tHolderList.add(tHolder);
			}

		}
		else if (isGoodBase(tTestBaseChar))
		{
			int iNewReadPos = iTrackIt + 1;

			this.tRegression.calcRowVector(iArrayOfQualityScores[iTrackIt], dAverageQualityScore, tTestBaseChar, iNewReadPos, 0, bNegativeStrand, iReadLength);

			dNewPhredScore = this.tRegression.runRegression();

			tHolder = RegressionScoreHolder.makeHolder(Character.toString(tTestBaseChar), new Float((float)dNewPhredScore), sBaseString, bIsBadBase, null);

			tHolderList.add(tHolder);
		}
		else
		{
			bIsBadBase = true;

			dNewPhredScore = iArrayOfQualityScores[iTrackIt];

			tHolder = RegressionScoreHolder.makeHolder(Character.toString(sBaseString.charAt(iTrackIt)), new Float((float)dOldPhredScore), sBaseString, bIsBadBase, Integer.valueOf(iArrayOfQualityScores[iTrackIt]));

			tHolderList.add(tHolder);
		}
			}

	public boolean doRegressionForBreak(String sBaseString, SortedBase tSortedBaseCurrent, int iTrackIt, char tTestBaseChar, int[] iArrayOfQualityScores, double dAverageQualityScore, boolean bNegativeStrand, int iReadLength, Iterator<SortedBase> tSortedBaseIterator, double dNewPhredScore, double dOldPhredScore, RegressionScoreHolder tHolder, List<RegressionScoreHolder> tHolderList, boolean bDoBreak)
			throws Exception
			{
		boolean bIsBadBase = false;

		tSortedBaseCurrent.add(iTrackIt, tTestBaseChar, iArrayOfQualityScores[iTrackIt], dAverageQualityScore, bNegativeStrand, iReadLength);
		if (tSortedBaseIterator.hasNext()) {
			tSortedBaseCurrent = (SortedBase)tSortedBaseIterator.next();
		}
		else if (isGoodBase(tTestBaseChar))
		{
			int iNewReadPos = iTrackIt + 1;

			this.tRegression.calcRowVector(iArrayOfQualityScores[iTrackIt], dAverageQualityScore, tTestBaseChar, iNewReadPos, 0, bNegativeStrand, iReadLength);

			dNewPhredScore = this.tRegression.runRegression();

			tHolder = RegressionScoreHolder.makeHolder(Character.toString(tTestBaseChar), new Float((float)dNewPhredScore), sBaseString, bIsBadBase, null);

			tHolderList.add(tHolder);
		}
		else
		{
			bIsBadBase = true;

			dNewPhredScore = iArrayOfQualityScores[iTrackIt];

			tHolder = RegressionScoreHolder.makeHolder(Character.toString(sBaseString.charAt(iTrackIt)), new Float((float)dOldPhredScore), sBaseString, bIsBadBase, new Integer(iArrayOfQualityScores[iTrackIt]));

			tHolderList.add(tHolder);
		}

		if (isGoodBase(tTestBaseChar))
		{
			int iNewReadPos = iTrackIt + 1;

			this.tRegression.calcRowVector(iArrayOfQualityScores[iTrackIt], dAverageQualityScore, tTestBaseChar, iNewReadPos, 0, bNegativeStrand, iReadLength);

			dNewPhredScore = this.tRegression.runRegression();

			tHolder = RegressionScoreHolder.makeHolder(Character.toString(tTestBaseChar), new Float((float)dNewPhredScore), sBaseString, bIsBadBase, null);

			tHolderList.add(tHolder);
		}
		else
		{
			bIsBadBase = true;

			dNewPhredScore = iArrayOfQualityScores[iTrackIt];

			tHolder = RegressionScoreHolder.makeHolder(Character.toString(sBaseString.charAt(iTrackIt)), new Float((float)dOldPhredScore), sBaseString, bIsBadBase, new Integer(iArrayOfQualityScores[iTrackIt]));

			tHolderList.add(tHolder);
		}

		return bDoBreak;
			}

	public void addRead(SAMRecord tSAMRecord, SAMFileWriter tOutputSAMWriter, SAMFileHeader tHeaderIn)
			throws IOException, Exception
			{
		this.tRecalScoresBuffer.setLength(0);
		this.tRecalScoresList.clear();

		List tHolderList = new ArrayList();

		String sChromosomeName = tSAMRecord.getReferenceName();
		int iReferenceLocation = tSAMRecord.getAlignmentStart();
		int iReadLength = tSAMRecord.getReadLength();

		shiftList(tSAMRecord, sChromosomeName, iReferenceLocation);

		byte[] bases = tSAMRecord.getReadBases();
		String basestring = StringUtil.bytesToString(bases);
		byte[] quals = tSAMRecord.getBaseQualities();

		int[] iquals = new int[quals.length];
		double avgq = 0.0D;
		for (int i = 0; i < quals.length; i++) {
			iquals[i] = quals[i];
			avgq += 1.0D * iquals[i];
		}

		avgq /= iquals.length;
		if (this.bApplyPhredScale) {
			avgq -= this.PHRED_SCALE_VALUE;
		}

		boolean negstrand = tSAMRecord.getReadNegativeStrandFlag();

		List abs = tSAMRecord.getAlignmentBlocks();

		Iterator tBlockIter = abs.iterator();

		Iterator iter = this.localregion.iterator();
		SortedBase current = (SortedBase)iter.next();
		int iAbsCounter = 0;

		boolean bFinishedReadString = false;

		int iTrackIt = 0;
		char tTestBaseChar;
		while (tBlockIter.hasNext())
		{
			AlignmentBlock ab = (AlignmentBlock)tBlockIter.next();

			iAbsCounter++;

			int rstart = ab.getReferenceStart();

			while (current.getRefLocation() < rstart) {
				current = (SortedBase)iter.next();
			}

			RegressionScoreHolder tHolder = null;

			double dNewPhredScore = 0.0D;
			double dOldPhredScore = 0.0D;

			bFinishedReadString = false;

			while (!bFinishedReadString)
			{
				if (basestring.length() - 1 <= iTrackIt) {
					bFinishedReadString = true;
				}
				else {
					for (int iTracker = 1; iTracker <= ab.getLength(); iTracker++)
					{
						if ((basestring.length() - 1 > iTrackIt) && (iTracker == ab.getLength()) && (!tBlockIter.hasNext()))
						{
							while (iTrackIt <= basestring.length() - 1)
							{
								if (tHolderList.size() == basestring.length())
								{
									break;
								}
								tTestBaseChar = basestring.charAt(iTrackIt);

								doRegression(basestring, current, iTrackIt, tTestBaseChar, iquals, avgq, negstrand, iReadLength, iter, dNewPhredScore, dOldPhredScore, tHolder, tHolderList);

								if (iTrackIt == basestring.length() - 1) {
									break;
								}
								iTrackIt++;
							}

							bFinishedReadString = true;
							break;
						}
						if (basestring.length() - 1 >= iTrackIt) {
							tTestBaseChar = basestring.charAt(iTrackIt);
							doRegressionForBreak(basestring, current, iTrackIt, tTestBaseChar, iquals, avgq, negstrand, iReadLength, iter, dNewPhredScore, dOldPhredScore, tHolder, tHolderList, true);

							iTrackIt++;

							if (tHolderList.size() == basestring.length())
							{
								break;
							}
							if ((iTracker == ab.getLength()) && (basestring.length() - 1 == iTracker) && (iTrackIt == iTracker)) {
								tTestBaseChar = basestring.charAt(iTrackIt);
								doRegressionForBreak(basestring, current, iTrackIt, tTestBaseChar, iquals, avgq, negstrand, iReadLength, iter, dNewPhredScore, dOldPhredScore, tHolder, tHolderList, true);

								break;
							}

							if ((iTracker == ab.getLength()) && (iAbsCounter > 1) && (basestring.length() - 1 == iTrackIt)) {
								tTestBaseChar = basestring.charAt(iTrackIt);
								doRegressionForBreak(basestring, current, iTrackIt, tTestBaseChar, iquals, avgq, negstrand, iReadLength, iter, dNewPhredScore, dOldPhredScore, tHolder, tHolderList, true);

								break;
							}

						}
						else if (basestring.length() - 1 <= iTrackIt) {
							if (tHolderList.size() == basestring.length()) {
								break;
							}
							tTestBaseChar = basestring.charAt(iTrackIt);
							bFinishedReadString = true;
							doRegressionForBreak(basestring, current, iTrackIt, tTestBaseChar, iquals, avgq, negstrand, iReadLength, iter, dNewPhredScore, dOldPhredScore, tHolder, tHolderList, true);

							break;
						}

					}

					bFinishedReadString = true;
					break;
				}

			}

		}

		Iterator tIter = tHolderList.iterator();

		while (tIter.hasNext())
		{
			RegressionScoreHolder tTempHolder = (RegressionScoreHolder)tIter.next();

			Float tFloat = tTempHolder.getScore();

			int iPhredScore = convertToPhredScale(tFloat);

			if (tIter.hasNext())
			{
				if ((tTempHolder.getBase().equals("N")) && (tTempHolder.isBadBase())) {
					this.tRecalScoresList.add(tTempHolder.getOriginalQualityScore());
				}
				else {
					this.tRecalScoresList.add(new Integer(iPhredScore));
				}

			}
			else if ((tTempHolder.getBase().equals("N")) && (tTempHolder.isBadBase())) {
				this.tRecalScoresList.add(tTempHolder.getOriginalQualityScore());
			}
			else {
				this.tRecalScoresList.add(new Integer(iPhredScore));
			}

		}

		byte[] tByteArray = new byte[this.tRecalScoresList.size()];
		Integer[] tIntegerArray = (Integer[])Arrays.copyOf(this.tRecalScoresList.toArray(), this.tRecalScoresList.toArray().length, Integer[].class);
		int iPos = 0;

		int d1 = tIntegerArray.length; 
		for (int dOldPhredScore = 0; dOldPhredScore < d1; dOldPhredScore++) { 
			Integer iScore =  tIntegerArray[dOldPhredScore];
			tByteArray[(iPos++)] = iScore.byteValue();
		}

		tSAMRecord.setBaseQualities(tByteArray);

		tOutputSAMWriter.addAlignment(tSAMRecord);

		this.tRecalScoresBuffer.setLength(0);
		this.tRecalScoresList.clear();
		tHolderList.clear();
			}

	private boolean isGoodBase(char charAt)
	{
		boolean bIsGoodBase = false;

		if (!this.tSetOfBadBases.contains(String.valueOf(charAt))) {
			bIsGoodBase = true;
		}
		return bIsGoodBase;
	}

	public static float byteArrayToFloat(byte[] test)
	{
		int bits = 0;
		int i = 0;
		for (int shifter = 3; shifter >= 0; shifter--) {
			bits |= (test[i] & 0xFF) << shifter * 8;
			i++;
		}
		return Float.intBitsToFloat(bits);
	}

	public static float byteArrayToFloat2(Byte[] test)
	{
		int bits = 0;
		int i = 0;
		for (int shifter = 3; shifter >= 0; shifter--) {
			Byte tByte = test[i];
			if (tByte != null) {
				bits |= (tByte.intValue() & 0xFF) << shifter * 8;
				i++;
			}
		}
		return Float.intBitsToFloat(bits);
	}

	public static byte[] intToByteArray(int param)
	{
		byte[] result = new byte[4];
		for (int i = 0; i < 4; i++) {
			int offset = (result.length - 1 - i) * 8;
			result[i] = (byte)(param >>> offset & 0xFF);
		}
		return result;
	}

	public static <T> T[] concat(T[] first, T[] second) {
		Object[] result = Arrays.copyOf(first, first.length + second.length);
		System.arraycopy(second, 0, result, first.length, second.length);
		return (T[]) result;
	}

	public static byte[] concatb(byte[] first, byte[] second) {
		byte[] result = Arrays.copyOf(first, first.length + second.length);
		System.arraycopy(second, 0, result, first.length, second.length);
		return result;
	}

	public boolean needsPhredScaling(SAMFileReader tFileReaderIn)
			throws IOException
			{
		boolean bNeedsPhredScaling = false;
		int iRecordsRead = 0;
		int iMinQualityScore = 0;

		for (SAMRecord tLoopRec : tFileReaderIn)
		{
			String chr = tLoopRec.getReferenceName();
			int refloc = tLoopRec.getAlignmentStart();

			shiftList(tLoopRec, chr, refloc);

			byte[] quals = tLoopRec.getBaseQualities();

			int[] iquals = new int[quals.length];

			for (int i = 0; i < quals.length; i++) {
				iquals[i] = quals[i];
			}

			Arrays.sort(iquals);
			iMinQualityScore = iquals[0];

			if ((iRecordsRead == this.RECORDS_TO_READ) && (iMinQualityScore >= this.PHRED_SCALE_VALUE)) {
				bNeedsPhredScaling = true;
				break;
			}if ((iRecordsRead == this.RECORDS_TO_READ) && (iMinQualityScore < this.PHRED_SCALE_VALUE))
			{
				break;
			}
			iRecordsRead++;
		}

		return bNeedsPhredScaling;
			}

	public void setPhredScaling(boolean bPhredScaleIn)
	{
		this.bApplyPhredScale = bPhredScaleIn;
		this.tRegression.setPhredScale(bPhredScaleIn);
	}

	private void shiftList(SAMRecord rec, String chr, int refloc)
			throws IOException
			{
		flushBefore(chr, refloc);

		if (this.localregion.size() != 0)
		{
			assert (((SortedBase)this.localregion.getFirst()).getRefLocation() == refloc);
		}

		if (this.localregion.size() == 0)
		{
			this.localregion.add(getSortedBase(chr, refloc));
		}

		int ll = ((SortedBase)this.localregion.getLast()).getRefLocation();
		while (ll <= rec.getAlignmentEnd())
		{
			ll++;
			this.localregion.addLast(getSortedBase(chr, ll));
		}
			}

	private SortedBase getSortedBase(String chrom, int ref_location)
	{
		return new SortedBase(chrom, ref_location);
	}

	public void flushAll() throws IOException
	{
		while (!this.localregion.isEmpty())
		{
			SortedBase localSortedBase = (SortedBase)this.localregion.pop();
		}
	}

	public int convertToPhredScale(Float tFloatIn)
	{
		int iPhredToReturn = 0;
		iPhredToReturn = (int)(-10.0D * Math.log10(tFloatIn.doubleValue()));
		return iPhredToReturn;
	}

	private class SortedBase
	{
		String chromosome = null;
		int refloc;
		List<Integer> read_locs = new ArrayList();
		List<Character> read_base = new ArrayList();
		List<Integer> base_qual = new ArrayList();
		List<Double> avgqual = new ArrayList();
		List<Boolean> reverse = new ArrayList();

		List<Integer> tReadLengthList = new ArrayList();

		public SortedBase(String chr, int rl)
		{
			this.chromosome = chr;
			this.refloc = rl;
		}
		public String getChromosome() {
			return this.chromosome; } 
		public int getRefLocation() { return this.refloc; }

		public void add(int readpos, char charAt, int b, double avgq, boolean r, int iReadLengthIn) {
			this.read_locs.add(Integer.valueOf(readpos));
			this.read_base.add(Character.valueOf(charAt));
			this.base_qual.add(Integer.valueOf(b));
			this.avgqual.add(Double.valueOf(avgq));
			this.reverse.add(Boolean.valueOf(r));

			this.tReadLengthList.add(Integer.valueOf(iReadLengthIn));
		}
	}
}