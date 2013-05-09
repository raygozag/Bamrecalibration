/*    */ package org.renci.sequencing.util;
/*    */ 
/*    */ public class RegressionScoreHolder
/*    */ {
/* 11 */   private String sBase = null;
/*    */ 
/* 13 */   private Float tFloat = Float.valueOf(0.0F);
/*    */ 
/* 15 */   private String sReadString = null;
/*    */   private boolean bIsBadBase;
/*    */   private Integer tOriginalQualityScore;
/*    */ 
/*    */   private RegressionScoreHolder(String sBaseIn, Float tRegressionScoreIn, String sReadStringIn, boolean bIsBadBaseIn, Integer tOriginalQualityScoreIn)
/*    */   {
/* 28 */     this.sBase = sBaseIn;
/* 29 */     this.tFloat = tRegressionScoreIn;
/* 30 */     this.sReadString = sReadStringIn;
/* 31 */     this.bIsBadBase = bIsBadBaseIn;
/* 32 */     this.tOriginalQualityScore = tOriginalQualityScoreIn;
/*    */   }
/*    */ 
/*    */   public static RegressionScoreHolder makeHolder(String sBaseIn, Float tRegressionScoreIn, String sReadStringIn, boolean bIsBadBaseIn, Integer tOriginalQualityScoreIn)
/*    */   {
/* 41 */     return new RegressionScoreHolder(sBaseIn, tRegressionScoreIn, sReadStringIn, bIsBadBaseIn, tOriginalQualityScoreIn);
/*    */   }
/*    */ 
/*    */   public boolean isBadBase()
/*    */   {
/* 52 */     return this.bIsBadBase;
/*    */   }
/*    */ 
/*    */   public Integer getOriginalQualityScore()
/*    */   {
/* 62 */     return this.tOriginalQualityScore;
/*    */   }
/*    */ 
/*    */   public String getBase()
/*    */   {
/* 72 */     return this.sBase;
/*    */   }
/*    */ 
/*    */   public Float getScore()
/*    */   {
/* 82 */     return this.tFloat;
/*    */   }
/*    */ 
/*    */   public String getReadString() {
/* 86 */     return this.sReadString;
/*    */   }
/*    */ }

/* Location:           /Users/raygoza/SparkleShare/courses/bioinformaticsI/project/ReQON/inst/java/BAMRecalibrator.jar
 * Qualified Name:     org.renci.sequencing.util.RegressionScoreHolder
 * JD-Core Version:    0.6.1
 */