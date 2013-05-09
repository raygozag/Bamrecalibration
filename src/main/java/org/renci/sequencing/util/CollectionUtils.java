package org.renci.sequencing.util;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Collections;
import java.util.Iterator;
import java.util.Vector;

import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;

public class CollectionUtils<T> {

	public static Vector<Integer> pmax(int a, Vector<Integer> v, double nraf){
		
		Vector<Integer> pmax= new Vector<Integer>();
		
		
		for(Integer cov: v) {
			
			if((double)a < (double)cov*nraf) {
				pmax.add(Integer.valueOf(a));
			}else {
				pmax.add(cov);
			}
			
			
		}
		
		
		
		return pmax;
	}
	
	
public static Vector<Double> pMax(int a, Vector<MutableInt> v, double nraf){
		
		Vector<Double> pmax= new Vector<Double>();
		double value= (double)a;
		
		for(MutableInt cov: v) {
			
			if(value > cov.doubleValue()*nraf) {
				pmax.add(value);
			}else {
				pmax.add(cov.doubleValue()*nraf);
			}
			
			
		}
		
		return pmax;
	}


public static Rle rle(Vector<Integer> v){
	
	Rle rle = new Rle();
	
	for(Integer i: v) {
		rle.Increase(i);
	}
	return rle;
}



public static Vector<Integer> filterBases(Vector<Integer> base_quality, Vector<Integer> pos){
	
	Vector<Integer> filtered_bases= new Vector<Integer>();
	
	for(int i=0; i < base_quality.size();i++) {
		if(base_quality.get(i)==1) {
			filtered_bases.add(pos.elementAt(i));
		}
	}
	
	
	return filtered_bases;
}


public static int sum(Vector<Integer> v) {
	int s=0;
	
	for(Integer inte:v) {
		s+=inte;
	}
	
	return s;
	
}


public static void deletePositions(Vector<Integer> d2,Vector<MutableInt> d3, Vector<String> d4, Vector<MutableInt> d5, Vector<MutableDouble> d6,Vector<String> d7,Vector<Integer> delete) {
	for(int i=0; i< delete.size(); i++) {
		int k =d2.indexOf(delete.elementAt(i));
		while(k>=0) {
		//	count.remove(k);
			//unique.remove(k);
			//d1.remove(k);
			d2.remove(k);
			d3.remove(k);
			d4.remove(k);
			d5.remove(k);
			d6.remove(k);
			//d7.remove(k);
			k = d2.indexOf(delete.elementAt(i),k);
		}
		
		
	}
}


public  Vector<T> removeAllIndexes(Vector<T> coll, Vector<Integer> idxs){
	
	Vector<T> integ = new Vector<T>();
	Vector<Integer> pos = new Vector<Integer>();
	pos.addAll(idxs);
	Collections.sort(pos);
	
	for(int i=0; i < pos.size();i++) {
		coll.setElementAt(null, pos.elementAt(i));
	}
	
	Iterator<T> it = coll.iterator();
	
	while(it.hasNext()) {
		if(it.next()==null) it.remove();
	}
	
	return coll;
}





public static void freeArray(Integer[] array) {
	for(int i=0; i < array.length; i++) {
		array[i]=null;
	}
}

public static int sum(Integer[] v) {
	int s=0;
	
	for(Integer inte:v) {
		s+=inte;
	}
	
	return s;
	
}


public static Rle filterPositions(Vector<Integer> Y, Vector<MutableInt> d3 ){
	Vector<Integer> filter = new Vector<Integer>();
	for(int i=0; i < Y.size(); i++) {
		if(Y.elementAt(i)==1) {
			filter.add(d3.elementAt(i).intValue());
		}
	}
	return rle(filter);
}

public static void writeIntVector(Vector<Integer> v,String filename) throws Exception {
	BufferedWriter wr = new BufferedWriter(new FileWriter(filename));
	
	for(Integer in: v) {
		wr.write(in.intValue()+"\n");
	}
	
	wr.close();
}


public static void writeMutableIntVector(Vector<MutableInt> v,String filename) throws Exception {
	BufferedWriter wr = new BufferedWriter(new FileWriter(filename));
	
	for(MutableInt in: v) {
		wr.write(in.intValue()+"\n");
	}
	
	wr.close();
}

public static void writeDoubleVector(Vector<Double> v,String filename) throws Exception {
	BufferedWriter wr = new BufferedWriter(new FileWriter(filename));
	
	for(Double in: v) {
		wr.write(in.doubleValue()+"\n");
	}
	
	wr.close();
}

public static void writeIntArray(Integer[] v,String filename) throws Exception {
	BufferedWriter wr = new BufferedWriter(new FileWriter(filename));
	
	for(Integer in: v) {
		wr.write(in.intValue()+"\n");
	}
	
	wr.close();
}

public static void writerle(Vector<Integer> pos, Vector<MutableInt> counts) throws Exception {
	
	BufferedWriter wr = new BufferedWriter(new FileWriter("countsrle.txt"));
	
	for(int i=0; i < pos.size(); i++) {
		wr.write(pos.elementAt(i)+"\t"+counts.elementAt(i)+"\n");
	}
	
	wr.close();
	
}

	
}
