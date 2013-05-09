/*     */ package org.renci.sequencing.util;
/*     */ 
/*     */ import java.io.BufferedWriter;
/*     */ import java.io.IOException;
/*     */ import java.io.PrintStream;
/*     */ import java.util.ArrayList;
/*     */ import java.util.List;
/*     */ 
/*     */ public class LogisticRegression
/*     */   implements IRegression
/*     */ {
/*     */   private List<String> tRowVectorList;
/*     */   private List<Float> tCoeffList;
/*     */   private List<Integer> tFlaggedList;
/*     */   private List<String> tReadFlagVectorList;
/*  25 */   private boolean bScalePhredBy33 = false;
/*     */ 
/*  27 */   private int PHRED_SCALE_VALUE = 33;
/*     */   private BufferedWriter tDebugWriter;
/*     */ 
/*     */   private LogisticRegression(List<Float> tCoeffListIn, List<Integer> tFlaggedListIn, boolean bIsPhredScaleBy33)
/*     */   {
/*  39 */     this.tCoeffList = tCoeffListIn;
/*  40 */     this.tFlaggedList = tFlaggedListIn;
/*  41 */     this.tReadFlagVectorList = new ArrayList();
/*  42 */     this.tRowVectorList = new ArrayList();
/*  43 */     this.bScalePhredBy33 = bIsPhredScaleBy33;
/*     */   }
/*     */ 
/*     */   private LogisticRegression(List<Float> tCoeffListIn, List<Integer> tFlaggedListIn, boolean bIsPhredScaleBy33, BufferedWriter tDebugWriterIn) {
/*  47 */     this.tCoeffList = tCoeffListIn;
/*  48 */     this.tFlaggedList = tFlaggedListIn;
/*  49 */     this.tReadFlagVectorList = new ArrayList();
/*  50 */     this.tRowVectorList = new ArrayList();
/*  51 */     this.bScalePhredBy33 = bIsPhredScaleBy33;
/*  52 */     this.tDebugWriter = tDebugWriterIn;
/*     */   }
/*     */ 
/*     */   public static LogisticRegression getInstance(List<Float> tCoeffList, List<Integer> tFlaggedList, boolean bIsPhredScaleBy33)
/*     */   {
/*  66 */     return new LogisticRegression(tCoeffList, tFlaggedList, bIsPhredScaleBy33);
/*     */   }
/*     */ 
/*     */   public static LogisticRegression getInstance(List<Float> tCoeffList, List<Integer> tFlaggedList, boolean bIsPhredScaleBy33, BufferedWriter tDebugWriterIn) {
/*  70 */     return new LogisticRegression(tCoeffList, tFlaggedList, bIsPhredScaleBy33, tDebugWriterIn);
/*     */   }
/*     */ 
/*     */   public double runRegression()
/*     */     throws Exception
/*     */   {
/*  76 */     double dZ = calculateZ(2, this.tCoeffList);
/*  77 */     double dResult = calculateFZ(dZ);
/*  78 */     return dResult;
/*     */   }
/*     */ 
/*     */   public void calcRowVector(int iPhredIn, double iAvgQCIn, char sBaseIn, int iReadPosIn, int iReadFlagIn, boolean bNegativeStrandIn, int iReadLengthIn)
/*     */   {
/* 101 */     this.tRowVectorList.clear();
/*     */ 
/* 103 */     this.tRowVectorList.add("1");
/*     */ 
/* 105 */     int iLocalPhred = 0;
/* 106 */     if (this.bScalePhredBy33)
/* 107 */       iLocalPhred = iPhredIn - this.PHRED_SCALE_VALUE;
/*     */     else {
/* 109 */       iLocalPhred = iPhredIn;
/*     */     }
/*     */ 
/* 112 */     this.tRowVectorList.add(Integer.toString(iLocalPhred));
/* 113 */     int iPhredIsZero = iPhredIn == 0 ? 1 : 0;
/* 114 */     this.tRowVectorList.add(Integer.toString(iPhredIsZero));
/* 115 */     this.tRowVectorList.add(Double.toString(iAvgQCIn));
/* 116 */     this.tRowVectorList.add(Integer.toString(getBaseValue(String.valueOf(sBaseIn), "A")));
/* 117 */     this.tRowVectorList.add(Integer.toString(getBaseValue(String.valueOf(sBaseIn), "C")));
/* 118 */     this.tRowVectorList.add(Integer.toString(getBaseValue(String.valueOf(sBaseIn), "G")));
/*     */ 
/* 121 */     if (bNegativeStrandIn)
/*     */     {
/* 123 */       int iReadPos = calcReadPos(iReadLengthIn, iReadPosIn);
/*     */ 
/* 125 */       if (this.tDebugWriter != null) {
/*     */         try {
/* 127 */           this.tDebugWriter.write(iPhredIn + "\t" + iAvgQCIn + "\t" + sBaseIn + "\t" + iReadPos + System.getProperty("line.separator"));
/*     */         } catch (IOException e) {
/* 129 */           e.printStackTrace();
/*     */         }
/*     */ 
/*     */       }
/*     */ 
/* 134 */       this.tRowVectorList.add(Integer.toString(iReadPos));
/* 135 */       this.tRowVectorList.addAll(getReadFlagVectorList(iReadPos));
/*     */     }
/*     */     else
/*     */     {
/* 139 */       if (this.tDebugWriter != null) {
/*     */         try {
/* 141 */           this.tDebugWriter.write(iPhredIn + "\t" + iAvgQCIn + "\t" + sBaseIn + "\t" + iReadPosIn + System.getProperty("line.separator"));
/*     */         } catch (IOException e) {
/* 143 */           e.printStackTrace();
/*     */         }
/*     */       }
/*     */ 
/* 147 */       this.tRowVectorList.add(Integer.toString(iReadPosIn));
/* 148 */       this.tRowVectorList.addAll(getReadFlagVectorList(iReadPosIn));
/*     */     }
/*     */   }
/*     */ 
/*     */   private List<String> getReadFlagVectorList(int iReadPosIn)
/*     */   {
/* 171 */     this.tReadFlagVectorList.clear();
/* 172 */     int iIndexOfFlaggedReadPos = this.tFlaggedList.indexOf(Integer.valueOf(iReadPosIn));
/* 173 */     if (iIndexOfFlaggedReadPos != -1)
/*     */     {
/* 175 */       for (int ii = 0; ii < this.tFlaggedList.size(); ii++) {
/* 176 */         if (ii == iIndexOfFlaggedReadPos)
/* 177 */           this.tReadFlagVectorList.add(String.valueOf(1));
/*     */         else
/* 179 */           this.tReadFlagVectorList.add(String.valueOf(0));
/*     */       }
/*     */     }
/*     */     else {
/* 183 */       for (int ii = 0; ii < this.tFlaggedList.size(); ii++) {
/* 184 */         this.tReadFlagVectorList.add(String.valueOf(0));
/*     */       }
/*     */     }
/*     */ 
/* 188 */     return this.tReadFlagVectorList;
/*     */   }
/*     */ 
/*     */   private int calcReadPos(int iReadLengthIn, int iReadPosIn)
/*     */   {
/* 202 */     int iReadPosCalc = 0;
/*     */ 
/* 204 */     iReadLengthIn++;
/* 205 */     iReadPosIn++;
/*     */ 
/* 209 */     iReadPosCalc = iReadLengthIn - iReadPosIn + 1;
/*     */ 
/* 211 */     return iReadPosCalc;
/*     */   }
/*     */ 
/*     */   private int getBaseValue(String sBaseValueIn, String sCompareStringIn)
/*     */   {
/* 242 */     int iReturnValue = 0;
/*     */ 
/* 244 */     iReturnValue = sBaseValueIn.equalsIgnoreCase(sCompareStringIn) ? 1 : 0;
/*     */ 
/* 246 */     return iReturnValue;
/*     */   }
/*     */ 
/*     */   private double calculateZ(int iBetaIntercept, List<Float> tCoeffListIn)
/*     */     throws Exception
/*     */   {
/* 256 */     double dZValue = 0.0D;
/* 257 */     dZValue = calcDotProduct(tCoeffListIn, this.tRowVectorList);
/* 258 */     return dZValue;
/*     */   }
/*     */ 
/*     */   private double calculateFZ(double iZValueIn)
/*     */   {
/* 268 */     double iDoubleReturn = 0.0D;
/*     */ 
/* 270 */     iDoubleReturn = 1.0D / (1.0D + Math.pow(2.718281828459045D, -iZValueIn));
/*     */ 
/* 272 */     return iDoubleReturn;
/*     */   }
/*     */ 
/*     */   private double calcDotProduct(List<Float> tCoeffListIn, List<String> tGeneratedListIn)
/*     */     throws Exception
/*     */   {
/* 281 */     double dZValue = 0.0D;
/* 282 */     if (tCoeffListIn.size() != tGeneratedListIn.size()) {
/* 283 */       printError("LogisticRegression::calcDotProduct: Input arrays not of same length: tCoeffListIn: " + tCoeffListIn.toString() + " tGeneratedListIn: " + tGeneratedListIn.toString());
/* 284 */       throw new Exception("LogisticRegression::calcDotProduct: Input arrays not of same length: tCoeffListIn: " + tCoeffListIn.toString() + " tGeneratedListIn: " + tGeneratedListIn.toString());
/*     */     }
/*     */ 
/* 288 */     for (int ii = 0; ii < tGeneratedListIn.size(); ii++) {
/* 289 */       dZValue += Double.valueOf((String)tGeneratedListIn.get(ii)).doubleValue() * ((Float)tCoeffListIn.get(ii)).floatValue();
/*     */     }
/*     */ 
/* 294 */     return dZValue;
/*     */   }
/*     */ 
/*     */   public void setPhredScale(boolean bPhredScaleIn)
/*     */   {
/* 304 */     this.bScalePhredBy33 = bPhredScaleIn;
/*     */   }
/*     */ 
/*     */   private void printError(String sErrorMessageIn)
/*     */   {
/* 314 */     System.out.println(sErrorMessageIn);
/*     */   }
/*     */ }

/* Location:           /Users/raygoza/SparkleShare/courses/bioinformaticsI/project/ReQON/inst/java/BAMRecalibrator.jar
 * Qualified Name:     org.renci.sequencing.util.LogisticRegression
 * JD-Core Version:    0.6.1
 */