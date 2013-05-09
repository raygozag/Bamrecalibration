package org.raygoza.sequencing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;
import java.util.Vector;
import java.util.zip.GZIPInputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.BooleanUtils;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.renci.sequencing.util.CollectionUtils;
import org.renci.sequencing.util.Rle;

import weka.classifiers.functions.Logistic;
import weka.core.Instances;





public class BamRecalibratorTrainer {

	public static void main(String[] args)  throws Exception{
		
		int nerr=2;
		double nraf = 0.05;
		
		Options cli_options = initOptions();
		
		CommandLineParser parser = new PosixParser();
		
		CommandLine cl = parser.parse(cli_options, args);
		
		List<String> argss = cl.getArgList();
		
		for(String k: argss) {
			System.out.println(k);
		}
		
		String data_file=cl.getOptionValue("read_data");
		if(data_file==null) {
			throw new Exception("The \"read_data\" parameter is required. Please specify it");
		}
		
		if(cl.hasOption("nerr")) {
			nerr = Integer.parseInt(cl.getOptionValue("nerr"));
		}
		
		
		Vector<String> d1 = new Vector<String>();
		Vector<Integer> d2 = new Vector<Integer>();
		Vector<MutableInt> d3 = new Vector<MutableInt>();
		Vector<String> d4 = new Vector<String>();
		Vector<MutableInt> d5 = new Vector<MutableInt>();
		Vector<MutableDouble> d6 = new Vector<MutableDouble>();
		Vector<String> d7 = new Vector<String>();
		Vector<Integer> pos= new Vector<Integer>();
		Vector<MutableInt> rle_unique_pos= new Vector<MutableInt>();
		Vector<Integer> matching = new Vector<Integer>();
		Vector<String> bases= new Vector<String>();
		bases.add("A");
		bases.add("C");
		bases.add("G");
		bases.add("T");
		
		
		
		
		System.out.println("Loading training data...");
		
		loadData(data_file,d1,d2,d3,d4,d5,d6,d7,pos,rle_unique_pos);
		System.out.println("Done loading training data...");
		processQualities(d5,d6);
		
		MutableInt startPos = Collections.min(d3);
		MutableInt endPos = Collections.max(d3);
		
		Integer readLen = endPos.intValue() - startPos.intValue() +1;
		
		processStartPosition(readLen,startPos.intValue(),d3,d7);
		
		
		
		filterdata(d1,d2,d3,d4,d5,d6,d7,pos,rle_unique_pos);
		
		TreeSet<Integer> uniquepos= new TreeSet<Integer>(d2);
		
		
		 
		 
		 System.out.println("Detecting sequencing errors");
		 
		matchUniqueToData(uniquepos,d2,matching);
		matching.add(Integer.valueOf(d2.size()));
		
		Vector<Integer> Y = new Vector<Integer>();
		Integer[] Useless = new Integer[d2.size()];
		for(int i=0; i < uniquepos.size(); i++ ) {
			
			Vector<Integer>  counts = getNucleotideCounts(matching.elementAt(i), matching.elementAt(i+1)-1, d4);
			
			
			int k= Collections.max(counts);
			
			
			for(int l=matching.elementAt(i); l <=  (matching.elementAt(i+1)-1);l++ ) {
				
				if(d4.elementAt(l).toUpperCase().equals("A") && counts.elementAt(0)==k) {
					Y.add(Integer.valueOf(0));
				}else if(d4.elementAt(l).toUpperCase().equals("C") && counts.elementAt(1)==k) {
					Y.add(Integer.valueOf(0));
				}else if(d4.elementAt(l).toUpperCase().equals("G") && counts.elementAt(2)==k) {
					Y.add(Integer.valueOf(0));
				}else if(d4.elementAt(l).toUpperCase().equals("T") && counts.elementAt(3)==k) {
					Y.add(Integer.valueOf(0));
				}else {
					Y.add(Integer.valueOf(1));
				}
				if(bases.contains(d4.elementAt(l).toUpperCase())) {
					Useless[l]= Integer.valueOf(0);
				}else {
					Useless[l]= Integer.valueOf(1);
				}
				
				
			}
			
			
		}
		
		
		
		System.out.println("Removing ambiguous positions");
		
		Rle coverage = CollectionUtils.rle(d2);
		
		Vector<Double> thresh = CollectionUtils.pMax(nerr, coverage.getCounts(), nraf);
		
		Vector<Integer> filtered = CollectionUtils.filterBases(Y, d2);
		
		Rle r = CollectionUtils.rle(filtered);
		
		Vector<Integer> amb =computeAmbiguous(coverage,r,thresh);
		
		Vector<Integer> amb_idx = filterAmbiguous(amb,Useless,d2);
		
		CollectionUtils<Integer> intutils= new CollectionUtils<Integer>();
		CollectionUtils<String>  strutils = new CollectionUtils<String>();
		CollectionUtils<MutableDouble> mdblUtils = new CollectionUtils<MutableDouble>();
		CollectionUtils<MutableInt> mintUtils = new CollectionUtils<MutableInt>();
		
	//	d1 = strutils.removeAllIndexes(d1, amb_idx);
		//d2 = intutils.removeAllIndexes(d2, amb_idx);
		
		d3 = mintUtils.removeAllIndexes(d3, amb_idx);
		
		d4 = strutils.removeAllIndexes(d4, amb_idx);
		d5 = mintUtils.removeAllIndexes(d5, amb_idx);
		d6 = mdblUtils.removeAllIndexes(d6, amb_idx);
		//d7 = strutils.removeAllIndexes(d7, amb_idx);
		d2.clear();
		d7.clear();
		Y = intutils.removeAllIndexes(Y, amb_idx);
		System.gc();
		//BufferedReader keyb= new BufferedReader(new InputStreamReader(System.in));
		//keyb.readLine();
		double err = (double)CollectionUtils.sum(Y)/Y.size();
		
		System.out.println("Training set error rate = "+err);
		System.out.println("Number of bases in filtered training set: "+Y.size());
		
		Rle hist =CollectionUtils.filterPositions(Y,d3);
		
		
		System.out.println("Finding flagged positions...");
		
		double constant = 1.5*((double) CollectionUtils.sum(Y)/(double)readLen);
		Vector<Integer> flagpos = new Vector<Integer>();
		for(int l=0; l < readLen; l++) {
			if(hist.getCounts().elementAt(l).doubleValue()> constant) {
				flagpos.add(hist.getElements().elementAt(l));
			}
		}
		Collections.sort(flagpos);
		
		int m = (int)Math.ceil(d3.size()/10000000.0);
		
		double[][] coeffs = new double[m][13];
		
		for(int j=0; j<m; j++) {
			Integer[] mark = new Integer[2];
			mark[0] = Math.round(((float)d3.size()/m)*(j));
			mark[1] = Math.min(Math.round((float)d3.size()/m)*(j+1), d3.size());
			//saveTemporaryDataset(mark[0],mark[1],d3,d4,d5,d6,Y,flagpos);
			Instances dataset = getWekaDataset(mark[0],mark[1],d3,d4,d5,d6,Y,flagpos);
			dataset.setClassIndex(dataset.numAttributes()-1);
			Logistic logistic = new Logistic();
			logistic.setRidge(0);
			logistic.buildClassifier(dataset);
			
			double[][] coeffx = logistic.coefficients();
			
			for(int i=0; i < coeffx.length; i++) {
				coeffs[j][i]=coeffx[i][0]*-1.0;
				System.out.println(coeffx[i][0]*-1.0);
			}
		}
		
		
		
		
		//writeIntArray(Y, "Y.txt");
		//writeDoubleVector(thresh, "thrs.txt");
		//writedata(d1,d2,d3,d4,d5,d6,d7);
		//writeIntVector(amb, "amb.txt");
		
		
		
	}
	
	
	public static Options initOptions() {
		Options opts = new Options();
		
		Option opt1 = new Option("read_data",true, "A tab delimited file with the bam file as text.");
		Option opt2 = new Option("ref_seq", "Path to reference genome.");
		Option opt3 = new Option("nerr", "maximum number of errors tolerated at a genomic position.");
		Option opt4 = new Option("nraf", "maximum non-reference allele frequency at a genomic position that is allowed.");
		Option opt5 = new Option("snp", "file of SNP locations to remove from training set before recalibration.");
		
		opts.addOption(opt1);
		opts.addOption(opt2);
		opts.addOption(opt3);
		opts.addOption(opt4);
		opts.addOption(opt5);
		
		return opts;
	}

	
	private static void filterdata( Vector<String> d1, Vector<Integer> d2,Vector<MutableInt> d3, Vector<String> d4, Vector<MutableInt> d5, Vector<MutableDouble> d6,Vector<String> d7,Vector<Integer> unique,Vector<MutableInt> count) {
		
		Vector<Integer> delete = new Vector<Integer>();
		
		for(int i=0; i < unique.size();i++) {
			if(count.elementAt(i).intValue()<3) {
				delete.add(unique.elementAt(i));
			}
		}
		
		CollectionUtils.deletePositions(d2, d3, d4, d5, d6, d7, delete);
		
		
	}
	
	private static void loadData(String fname, Vector<String> d1, Vector<Integer> d2,Vector<MutableInt> d3, Vector<String> d4, Vector<MutableInt> d5, Vector<MutableDouble> d6,Vector<String> d7,Vector<Integer> unique,Vector<MutableInt> count) throws Exception
	{
		
		BufferedReader rd = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fname))),20*1024*1024);
		
		String line="";
		Hashtable<String,String> bases= new Hashtable<String,String>();
		bases.put("A","A");
		bases.put("C","C");
		bases.put("G","G");
		bases.put("T","G");
		bases.put("N", "N");
		bases.put("+", "+");
		bases.put("-", "-");
		String chrom="";
		Vector<Integer> un = new Vector<Integer>();
		int i=-1;
		while(true) {
			line = rd.readLine();
			if(line==null) break;
			
			String[] vals = line.split("\t");
			
			if(!chrom.equals(vals[0].trim())) {
				chrom=vals[0];
			}
			
			//d1.add(vals[0]);
			
			Integer pos = Integer.valueOf(vals[1]);
			d2.add(pos);
			d3.add(new MutableInt(Integer.parseInt(vals[2])));
			d4.add(bases.get(vals[3].toUpperCase()));
			d5.add(new MutableInt(Integer.parseInt(vals[4])));
			d6.add(new MutableDouble(Double.parseDouble(vals[5])));
			d7.add(bases.get(vals[6]));
			

			if(!unique.contains(pos)) {
				unique.add(pos);
				count.addElement(new MutableInt(1));
				i++;
			}else {
				count.elementAt(i).increment();
			}
			
			
		}
		
		
		
	}

	
	private static void matchUniqueToData(TreeSet<Integer> unique, Vector<Integer> d2,Vector<Integer> match) {
		int from=0;
		for(Integer value: unique) {
			int idx = d2.indexOf(value	,from);
			from=idx;
			match.addElement(Integer.valueOf(idx));
		}
		
		
	}
	
	private static void processQualities(Vector<MutableInt> d5, Vector<MutableDouble> d6) {
		
		
		MutableInt mind5 = Collections.min(d5);
		
		if(mind5.intValue()>32) {
			for(int i=0; i < d5.size(); i++) {
				d5.elementAt(i).setValue(d5.elementAt(i).intValue()-33);
				d6.elementAt(i).setValue(d6.elementAt(i).doubleValue()-33.0);
			}
		}
		
		
		
	}
	
	
	private static void processStartPosition(int readlen, int startPos, Vector<MutableInt> d3, Vector<String> d7) {
		
		int add=0;
		
		if(startPos==0) {
			add=1;
		}
		
		for(int i=0; i< d3.size(); i++) {
		
			d3.elementAt(i).setValue(Math.abs((d3.elementAt(i).intValue()+add)-(readlen+1)*(d7.elementAt(i).equals("-")? 1:0)));
			
			
		}
		
		
		
	}
	
public static Instances getWekaDataset(int start, int end,Vector<MutableInt> d3, Vector<String> d4, Vector<MutableInt> d5,Vector<MutableDouble> d6,Vector<Integer> Y,Vector<Integer> flagged)  throws Exception{
		
		Random rnd = new Random();
		String filename = "/tmp/test"+rnd.nextInt()+".arff";
		StringWriter stwr= new StringWriter();
		BufferedWriter wr = new BufferedWriter(stwr);
		int isA=0;
		int isC=0;
		int isG=0;
		int[] flags = new int[flagged.size()];
		 DecimalFormat df = new DecimalFormat("#.######");
		 wr.write("@relation test\n" + 
		 		"\n" + 
		 		"@attribute a1 numeric\n" + 
		 		"@attribute a2 numeric\n" + 
		 		"@attribute a3 numeric\n" + 
		 		"@attribute a4 numeric\n" + 
		 		"@attribute a5 numeric\n" + 
		 		"@attribute a6 numeric\n" + 
		 		"@attribute a7 numeric\n" + 
		 		"@attribute a8 numeric\n" + 
		 		"@attribute a9 numeric\n" + 
		 		"@attribute a10 numeric\n" + 
		 		"@attribute a11 numeric\n" + 
		 		"@attribute a12 numeric\n" + 
		 		"@attribute a13 numeric\n" + 
		 		"@attribute a14 {0,1}\n" + 
		 		"\n" + 
		 		"@data\n");
		 
		for(int i=start; i <end; i++) {
			int d5e = d5.elementAt(i).intValue();
			if(d4.elementAt(i).equals("A")) {
				isA=1;
				isC=0;
				isG=0;
			}else if(d4.elementAt(i).equals("C")) {
				isA=0;
				isC=1;
				isG=0;
			}else if(d4.elementAt(i).equals("G")) {
				isA=0;
				isC=0;
				isG=1;
			}else {
				isA=0;
				isC=0;
				isG=0;
			}
			
			int flag = flagged.indexOf(d3.elementAt(i).intValue());
			
			
			for(int l=0; l < flags.length; l++) {
				flags[l]=0;
			}
			if(flag!=-1) {
				flags[flag]=1;
			}
			
			wr.write("1,"+d5e+ ","+BooleanUtils.toInteger(d5e==0)+","+df.format(d6.elementAt(i).doubleValue())+","+isA+","+isC+","+isG+","+d3.elementAt(i)+","+flags[0]);
			
			for(int j=1; j < flags.length;j++) {
				wr.write(","+flags[j]);
			}
			wr.write(","+Y.elementAt(i)+"\n");
			
		}
		
		wr.close();
		
		Instances data = new Instances(new StringReader(stwr.toString()));
		System.gc();
		return data;
		
	}
	

	private static void writedata(Vector<String> d1, Vector<Integer> d2,Vector<MutableInt> d3, Vector<String> d4, Vector<MutableInt> d5, Vector<MutableDouble> d6,Vector<String> d7) throws Exception{
		
		BufferedWriter wr = new BufferedWriter(new FileWriter("debug.txt"));
		
		wr.write(d2.elementAt(0)+"\t"+d3.elementAt(0)+"\t"+d4.elementAt(0)+"\t"+d5.elementAt(0)+"\t"+d6.elementAt(0)+"\t"+d7.elementAt(0)+"\n");
		
		
		for(int i=1; i < d1.size(); i++) {
			wr.write(d1.elementAt(i)+"\t"+d2.elementAt(i)+"\t"+d3.elementAt(i)+"\t"+d4.elementAt(i)+"\t"+d5.elementAt(i)+"\t"+d6.elementAt(i)+"\t"+d7.elementAt(i)+"\n");
		}
		
		wr.close();
	}
	
	private static Vector<Integer> getNucleotideCounts(int start, int end, Vector<String> sequence){
		
		int A=0;
		int C=0;
		int G=0;
		int T=0;
		
		for(int i=start;i <= end; i++) {
			if(sequence.elementAt(i).toUpperCase().equals("A")) {
				A++;
				continue;
			}
			if(sequence.elementAt(i).toUpperCase().equals("C")) {
				C++;
				continue;
			}
			if(sequence.elementAt(i).equals("G")) {
				G++;
				continue;
			}
			if(sequence.elementAt(i).equals("T")) {
				T++;
				continue;
			}
		}
		
		
		
		Vector<Integer> counts = new Vector<Integer>();
		counts.add(Integer.valueOf(A));
		counts.add(Integer.valueOf(C));
		counts.add(Integer.valueOf(G));
		counts.add(Integer.valueOf(T));
		
		return counts;
	}

	private static Vector<Integer> computeAmbiguous(Rle coverage, Rle r, Vector<Double> thresh){
		
		Vector<Integer> pos = new Vector<Integer>();
		
		
		Vector<Integer> elco = coverage.getElements();
		Vector<Integer> elr = r.getElements();
		Vector<MutableInt> cor = r.getCounts();
		
		
		for(int i=0; i < elco.size(); i++){
			
			if(elr.contains(elco.get(i))) {
				int k= elr.indexOf(elco.get(i));
				if(cor.get(k).doubleValue()>thresh.elementAt(i)) {
					pos.add(elr.elementAt(k));
				}
				
			}
		}
		
		
		return pos;
		
	}
	
	public static Vector<Integer> filterAmbiguous(Vector<Integer> amb,Integer[] Useless,Vector<Integer> d2) {
		
		Vector<Integer> pos = new Vector<Integer>();
		
		for(int i=0; i < d2.size(); i++) {
			if(Useless[i]==1) {
				pos.add(i);
				
				Useless[i]=null;
				continue;
			}
			
			if(amb.contains(d2.elementAt(i))) {
				
				pos.add(i);
			}	
			Useless[i]=null;
		}
		
		
		 
		return pos;
		
	}
	
	public static void saveTemporaryDataset(int start, int end,Vector<MutableInt> d3, Vector<String> d4, Vector<MutableInt> d5,Vector<MutableDouble> d6,Vector<Integer> Y,Vector<Integer> flagged)  throws Exception{
		
		BufferedWriter wr = new BufferedWriter(new FileWriter("/Volumes/Raygoza2/test.txt"));
		int isA=0;
		int isC=0;
		int isG=0;
		int[] flags = new int[flagged.size()];
		 DecimalFormat df = new DecimalFormat("#.######");
		 
		for(int i=start; i <end; i++) {
			int d5e = d5.elementAt(i).intValue();
			if(d4.elementAt(i).equals("A")) {
				isA=1;
				isC=0;
				isG=0;
			}else if(d4.elementAt(i).equals("C")) {
				isA=0;
				isC=1;
				isG=0;
			}else if(d4.elementAt(i).equals("G")) {
				isA=0;
				isC=0;
				isG=1;
			}else {
				isA=0;
				isC=0;
				isG=0;
			}
			
			int flag = flagged.indexOf(d3.elementAt(i).intValue());
			
			
			for(int l=0; l < flags.length; l++) {
				flags[l]=0;
			}
			if(flag!=-1) {
				flags[flag]=1;
			}
			
			wr.write("1\t"+d5e+ "\t"+BooleanUtils.toInteger(d5e==0)+"\t"+df.format(d6.elementAt(i).doubleValue())+"\t"+isA+"\t"+isC+"\t"+isG+"\t"+d3.elementAt(i)+"\t"+flags[0]);
			
			for(int j=1; j < flags.length;j++) {
				wr.write("\t"+flags[j]);
			}
			wr.write("\t"+Y.elementAt(i)+"\n");
			
		}
		
		wr.close();
		
	}
	
	
	
}
