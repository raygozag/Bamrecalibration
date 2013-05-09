/*     */ package org.renci.sequencing.util;
/*     */ 
/*     */ import java.io.BufferedReader;
/*     */ import java.io.BufferedWriter;
/*     */ import java.io.File;
/*     */ import java.io.FileNotFoundException;
/*     */ import java.io.FileReader;
/*     */ import java.io.FileWriter;
/*     */ import java.io.IOException;
/*     */ import java.io.PrintStream;
/*     */ import java.util.ArrayList;
/*     */ import java.util.Date;
/*     */ import java.util.HashSet;
/*     */ import java.util.Iterator;
/*     */ import java.util.List;
/*     */ import java.util.Set;
/*     */ import net.sf.picard.cmdline.CommandLineProgram;
/*     */ import net.sf.picard.cmdline.Option;
/*     */ import net.sf.picard.io.IoUtil;
/*     */ import net.sf.samtools.SAMFileHeader;
/*     */ import net.sf.samtools.SAMFileHeader.SortOrder;
/*     */ import net.sf.samtools.SAMFileReader;
/*     */ import net.sf.samtools.SAMFileReader.ValidationStringency;
/*     */ import net.sf.samtools.SAMFileWriter;
/*     */ import net.sf.samtools.SAMFileWriterFactory;
/*     */ import net.sf.samtools.SAMRecord;
/*     */ import net.sf.samtools.SAMRecordIterator;
/*     */ import net.sf.samtools.SAMSequenceDictionary;
/*     */ import net.sf.samtools.SAMSequenceRecord;
/*     */ 
/*     */ public class BAMRecalibrator extends CommandLineProgram
/*     */ {
/*     */ 
/*     */   @Option(shortName="I", doc="The input BAM file to process")
/*     */   public File INPUT;
/*     */ 
/*     */   @Option(doc="The name of the BAM output file.")
/*     */   public File BAM_OUTPUT;
/*     */ 
/*     */   @Option(doc="The name of the coefficients text file.")
/*     */   public File COEFF;
/*     */ 
/*     */   @Option(doc="The name of the flagged positions file.")
/*     */   public File FLAGGED;
/*     */ 
/*     */   @Option(doc="If set to true calculate mean quality over aligned reads only")
/* 113 */   public boolean ALIGNED_READS_ONLY = true;
/*     */ 
/*     */   @Option(shortName="PF", doc="If set to true calculate mean quality over PF reads only")
/* 116 */   public boolean PF_READS_ONLY = true;
/*     */ 
/*     */   @Option(doc="If set to true, include quality for no-call bases in the distribution")
/* 119 */   public boolean INCLUDE_NO_CALLS = true;
/*     */ 
/*     */   @Option(shortName="S", doc="If set to true, adjust the phred scale by 33")
/* 122 */   public boolean AUTO_PHRED_SCALE = true;
/*     */ 
/*     */   @Option(doc="Stop after processing N reads. Mostly for debugging purposes.")
/* 125 */   public int STOP_AFTER = 0;
/*     */ 
/*     */   @Option(doc="Region to unpack of the form 1:1-2500000 'chromosome:startpos-endpos'", optional=true)
/* 128 */   public String REGION = null;
/*     */ 
/*     */   @Option(doc="If set to true, display available chromome sequence names and sequence lengths for the input BAM file.", mutex={"BAM_OUTPUT", "COEFF", "FLAGGED"})
/* 131 */   public boolean TARGETS = false;
/*     */ 
/*     */   @Option(doc="Debug mode")
/* 134 */   public boolean DEBUG = false;
/*     */ 
/*     */   @Option(doc="ignore")
/* 137 */   public boolean IGNORE_WARNINGS = true;
/*     */ 
/* 140 */   private boolean bHasRegion = false;
/*     */   private int iStartPos;
/*     */   private int iEndPos;
/*     */   private String sChromosome;
/*     */   private Set<String> tSetOfBadBases;
/*     */   private static SAMFileWriter outputSam;
/*     */ 
/*     */   public BAMRecalibrator()
/*     */   {
/* 160 */     this.QUIET = Boolean.valueOf(true);
/*     */   }
/*     */ 
/*     */   public static void main(String[] args)
/*     */   {
/*     */     try
/*     */     {
/* 168 */       if (runIt(args) != 0) {
/* 169 */         System.out.println("Exiting with an error and a non-zero status.");
/* 170 */         throw new Exception("Exiting with an error and a non-zero status");
/*     */       }
/*     */ 
/*     */     }
/*     */     catch (RuntimeException rte)
/*     */     {
/* 177 */       outputSam.close();
/* 178 */       rte.printStackTrace();
/* 179 */       doRuntimeExceptionCIGARErrorAndExit(rte);
/*     */     }
/*     */     catch (Exception e) {
/* 182 */       System.err.println("Caught an exception : " + e.getMessage());
/* 183 */       e.printStackTrace();
/* 184 */       outputSam.close();
/*     */     }
/*     */   }
/*     */ 
/*     */   private static int runIt(String[] args)
/*     */   {
/* 204 */     int iReturn = new BAMRecalibrator().instanceMain(args);
/* 205 */     if (iReturn == 0) {
/* 206 */       System.out.println("[" + new Date() + "] " + "org.renci.sequencing.util.BAMRecalibrator done");
/* 207 */       System.out.println("Runtime.totalMemory()=" + Runtime.getRuntime().totalMemory());
/*     */     }
/* 209 */     return iReturn;
/*     */   }
/*     */ 
/*     */   public static void doRuntimeExceptionCIGARErrorAndExit(Exception tException)
/*     */   {
/* 225 */     if (tException.getMessage().contains("CIGAR")) {
/* 226 */       System.out.println("hey hey hey ###");
/* 227 */     } else if (tException.getMessage().contains("SAM validation error")) {
/* 228 */       tException.printStackTrace();
/* 229 */       System.err.println("CIGAR should have zero elements for unmapped string.  Maybe aligned using BWA?  Exiting recalibration.");
/*     */     }
/*     */     else {
/* 232 */       tException.printStackTrace();
/* 233 */       System.err.println("Caught a RuntimeException: " + tException.getMessage());
/*     */     }
/*     */   }
/*     */ 
/*     */   protected int doWork()
/*     */   {
/*     */     try
/*     */     {
/* 247 */       this.tSetOfBadBases = makeBadBaseSet();
/*     */ 
/* 250 */       SAMFileReader tSAMFileReaderIn = new SAMFileReader(this.INPUT);
/*     */ 
/* 252 */       SAMFileReader in2 = null;
/*     */ 
/* 254 */       tSAMFileReaderIn.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
/*     */ 
/* 257 */       if (this.TARGETS) {
/* 258 */         handleTargets(tSAMFileReaderIn);
/*     */       }
/*     */ 
/* 261 */       IoUtil.assertFileIsReadable(this.INPUT);
/*     */ 
/* 264 */       List tCoeffList = loadCoeffArrayFromFile(this.COEFF);
/* 265 */       List tFlaggedList = loadFlaggedArrayFromFile(this.FLAGGED);
/*     */ 
/* 268 */       SAMFileHeader tTestHeader = tSAMFileReaderIn.getFileHeader();
/* 269 */       tTestHeader.addComment("Quality scores were recalibrated with ReQON.");
/*     */ 
/* 271 */       outputSam = new SAMFileWriterFactory().makeBAMWriter(tTestHeader, true, this.BAM_OUTPUT);
/*     */ 
/* 274 */       if (this.REGION != null) {
/* 275 */         handleRegionParameters();
/*     */       }
/*     */ 
/* 278 */       if (!tSAMFileReaderIn.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
/* 279 */         System.out.println("WARNING: This BAM's header does not say it is sorted (perhaps it was sorted by samtools?)");
/*     */       }
/*     */ 
/* 283 */       SortedBaseHolder sbh = null;
/* 284 */       if (this.DEBUG) {
/* 285 */         String sDebugOutputDirectory = this.BAM_OUTPUT.getParent();
/* 286 */         String sDebugOutputFileAndPath = sDebugOutputDirectory + File.separator + "debug.txt";
/* 287 */         sbh = new SortedBaseHolder(tCoeffList, tFlaggedList, this.tSetOfBadBases, this.AUTO_PHRED_SCALE, new BufferedWriter(new FileWriter(sDebugOutputFileAndPath)));
/*     */       } else {
/* 289 */         sbh = new SortedBaseHolder(tCoeffList, tFlaggedList, this.tSetOfBadBases, this.AUTO_PHRED_SCALE);
/*     */       }
/*     */ 
/* 292 */       if ((this.AUTO_PHRED_SCALE) && (sbh.needsPhredScaling(tSAMFileReaderIn)))
/* 293 */         sbh.setPhredScaling(true);
/*     */       else {
/* 295 */         sbh.setPhredScaling(false);
/*     */       }
/* 297 */       tSAMFileReaderIn.close();
/*     */ 
/* 299 */       tSAMFileReaderIn = new SAMFileReader(this.INPUT);
/*     */ 
/* 302 */       if (this.bHasRegion)
/*     */       {
/* 304 */         if ((tSAMFileReaderIn == null) || (outputSam == null) || (sbh == null)) System.out.println("SAMFileReader is null or SAMOutputWriter is null or SortedBaseHolder is null.");
/* 305 */         processRegionQuery(tSAMFileReaderIn, outputSam, sbh);
/*     */       }
/*     */       else
/*     */       {
/* 310 */         int processed = 0;
/*     */ 
/* 314 */         for (SAMRecord rec : tSAMFileReaderIn)
/*     */         {
/* 316 */           if ((!this.PF_READS_ONLY) || (!rec.getReadFailsVendorQualityCheckFlag()))
/*     */           {
/* 319 */             sbh.addRead(rec, outputSam, outputSam.getFileHeader());
/* 320 */             if (this.STOP_AFTER > 0) { processed++; if (processed > this.STOP_AFTER) break; }
/*     */           }
/*     */         }
/* 323 */         outputSam.close();
/*     */       }
/*     */     }
/*     */     catch (IOException e) {
/* 327 */       System.out.println("I/O Error.  Underlying cause:\n " + e.getMessage());
/* 328 */       e.printStackTrace();
/*     */     }
/*     */     catch (RuntimeException e) {
/* 331 */       e.printStackTrace();
/* 332 */       doRuntimeExceptionCIGARErrorAndExit(e);
/*     */     }
/*     */     catch (Exception e) {
/* 335 */       e.printStackTrace();
/*     */     }
/*     */ 
/* 339 */     return 0;
/*     */   }
/*     */ 
/*     */   public void printAttribute(SAMFileReader tSAMFileReader, String sAttributeIn)
/*     */   {
/* 350 */     tSAMFileReader = new SAMFileReader(this.BAM_OUTPUT);
/* 351 */     SAMRecordIterator tI2 = tSAMFileReader.iterator();
/* 352 */     int iCounter = 0;
/* 353 */     while (tI2.hasNext()) {
/* 354 */       iCounter++;
/*     */ 
/* 356 */       SAMRecord tRec = (SAMRecord)tI2.next();
/*     */ 
/* 359 */       System.out.println("tRec.getAttribute(Q2): " + tRec.getAttribute("Q2"));
/* 360 */       System.out.println("tRec.getReadString() " + tRec.getReadString());
/* 361 */       System.out.println();
/*     */     }
/*     */   }
/*     */ 
/*     */   public Set<String> makeBadBaseSet()
/*     */   {
/* 374 */     Set tBadBaseSet = new HashSet();
/* 375 */     tBadBaseSet.clear();
/*     */ 
/* 377 */     tBadBaseSet.add("M");
/* 378 */     tBadBaseSet.add("R");
/* 379 */     tBadBaseSet.add("S");
/* 380 */     tBadBaseSet.add("V");
/* 381 */     tBadBaseSet.add("W");
/* 382 */     tBadBaseSet.add("Y");
/* 383 */     tBadBaseSet.add("H");
/* 384 */     tBadBaseSet.add("K");
/* 385 */     tBadBaseSet.add("D");
/* 386 */     tBadBaseSet.add("B");
/* 387 */     tBadBaseSet.add("N");
/*     */ 
/* 389 */     return tBadBaseSet;
/*     */   }
/*     */ 
/*     */   private void handleTargets(SAMFileReader tSAMFileReaderIn)
/*     */   {
/* 400 */     SAMFileHeader tHeader = tSAMFileReaderIn.getFileHeader();
/* 401 */     SAMSequenceDictionary tDict = tHeader.getSequenceDictionary();
/* 402 */     List tSRecordList = tDict.getSequences();
/* 403 */     Iterator tIter3 = tSRecordList.iterator();
/* 404 */     System.out.println("Sequence Name:\t Sequence Length");
/* 405 */     while (tIter3.hasNext()) {
/* 406 */       SAMSequenceRecord tSSRecord = (SAMSequenceRecord)tIter3.next();
/* 407 */       System.out.println(tSSRecord.getSequenceName() + ":\t " + tSSRecord.getSequenceLength());
/*     */     }
/*     */   }
/*     */ 
/*     */   private void handleRegionParameters()
/*     */     throws Exception
/*     */   {
/*     */     try
/*     */     {
/* 422 */       String[] sRegionArray = this.REGION.split(":");
/* 423 */       if (sRegionArray.length != 2) {
/* 424 */         System.out.println("REGION parameter is not of the form 1:1-2500000 where '1' is the chromosome, '1-2500000' is the start and end position.");
/* 425 */         throw new Exception("REGION parameter is not of the form 1:1-2500000 where '1' is the chromosome, '1-2500000' is the start and end position.");
/*     */       }
/*     */ 
/* 428 */       this.sChromosome = sRegionArray[0];
/* 429 */       String[] sPositionArray = sRegionArray[1].split("-");
/* 430 */       if (sPositionArray.length != 2) {
/* 431 */         System.out.println("REGION parameter does not include a start and end position.");
/* 432 */         throw new Exception("REGION parameter does not include a start and end position.");
/*     */       }
/*     */ 
/* 435 */       this.iStartPos = Integer.parseInt(sPositionArray[0]);
/* 436 */       this.iEndPos = Integer.parseInt(sPositionArray[1]);
/*     */ 
/* 440 */       this.bHasRegion = true;
/*     */     } catch (NumberFormatException e) {
/* 442 */       e.printStackTrace();
/*     */     }
/*     */   }
/*     */ 
/*     */   private void processRegionQuery(SAMFileReader tSAMFileReaderIn, SAMFileWriter outputSam, SortedBaseHolder sbh)
/*     */     throws Exception
/*     */   {
/*     */     try
/*     */     {
/* 455 */       int processed = 0;
/*     */ 
/* 457 */       SAMRecordIterator tSRIter = tSAMFileReaderIn.query(this.sChromosome, this.iStartPos, this.iEndPos, false);
/* 458 */       SAMFileHeader tHeader = new SAMFileHeader();
/* 459 */       while (tSRIter.hasNext()) {
/* 460 */         SAMRecord tRecord = (SAMRecord)tSRIter.next();
/*     */ 
/* 462 */         if ((!this.PF_READS_ONLY) || (!tRecord.getReadFailsVendorQualityCheckFlag()))
/*     */         {
/* 464 */           if ((tRecord == null) || (outputSam == null) || (sbh == null)) System.out.println("SAMRecord is null.");
/*     */ 
/* 467 */           sbh.addRead(tRecord, outputSam, tHeader);
/*     */ 
/* 469 */           if (this.STOP_AFTER > 0) { processed++; if (processed > this.STOP_AFTER) break; }
/*     */         }
/*     */       }
/* 472 */       outputSam.close();
/* 473 */       tSAMFileReaderIn.close();
/*     */     }
/*     */     catch (IOException e)
/*     */     {
/* 477 */       e.printStackTrace();
/*     */     }
/*     */   }
/*     */ 
/*     */   private List<Integer> loadFlaggedArrayFromFile(File tFileIn)
/*     */   {
/* 489 */     List tList = new ArrayList();
/*     */     try {
/* 491 */       FileReader tInStream = new FileReader(tFileIn);
/* 492 */       BufferedReader tBInStream = new BufferedReader(tInStream);
/* 493 */       String sLine = null;
/* 494 */       while ((sLine = tBInStream.readLine()) != null)
/* 495 */         tList.add(Integer.valueOf(sLine));
/*     */     }
/*     */     catch (FileNotFoundException e) {
/* 498 */       e.printStackTrace();
/*     */     }
/*     */     catch (IOException e) {
/* 501 */       e.printStackTrace();
/*     */     }
/*     */ 
/* 505 */     return tList;
/*     */   }
/*     */ 
/*     */   private List<Float> loadCoeffArrayFromFile(File tFileIn)
/*     */   {
/* 514 */     List tList = new ArrayList();
/*     */     try
/*     */     {
/* 517 */       FileReader tInStream = new FileReader(tFileIn);
/* 518 */       BufferedReader tBInStream = new BufferedReader(tInStream);
/* 519 */       String sLine = null;
/* 520 */       while ((sLine = tBInStream.readLine()) != null)
/* 521 */         tList.add(Float.valueOf(sLine));
/*     */     }
/*     */     catch (FileNotFoundException e) {
/* 524 */       e.printStackTrace();
/*     */     } catch (IOException e) {
/* 526 */       e.printStackTrace();
/*     */     }
/*     */ 
/* 529 */     return tList;
/*     */   }
/*     */ }

/* Location:           /Users/raygoza/SparkleShare/courses/bioinformaticsI/project/ReQON/inst/java/BAMRecalibrator.jar
 * Qualified Name:     org.renci.sequencing.util.BAMRecalibrator
 * JD-Core Version:    0.6.1
 */