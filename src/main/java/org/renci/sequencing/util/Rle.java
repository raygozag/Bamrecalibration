package org.renci.sequencing.util;

import java.util.Vector;

import org.apache.commons.lang3.mutable.MutableInt;

public class Rle {

	Vector<Integer> elements=new Vector<Integer>();
	Vector<MutableInt> counts= new Vector<MutableInt>();
	
	
	public void Rle() {
		
		
	}
	
	public void Increase(Integer i) {
		
		if(!elements.contains(i)) {
			elements.add(i);
			counts.add(new MutableInt(1));
			return;
		}
		
		int k = elements.indexOf(i);
		
		counts.elementAt(k).increment();
		
		
		
	}
	
	
	public int getCountFor(Integer i) {
		int k = elements.indexOf(i);
		
		return counts.elementAt(k).intValue();
		
	}

	public Vector<Integer> getElements() {
		return elements;
	}

	
	public Vector<MutableInt> getCounts() {
		return counts;
	}

	
	
}

