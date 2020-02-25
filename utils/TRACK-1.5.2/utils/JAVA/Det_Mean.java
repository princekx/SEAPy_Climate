package tigge;
import java.util.*;
import java.math.*;
import java.io.*;
import EPSAnalysis.*;

public class Det_Mean{

        private static final int NUM_TRACKS_MAX = 1000000;

	public static void main(String[] args)throws IOException{
		
		
		String date = null;
		String filename = null;
		
		String fptmp = null;
		
		String[] fpAnalysis = new String[NUM_TRACKS_MAX];
		
		double[] tmpInt1 = new double[Constants.LT_MAX];
		double[] sumInt1 = new double[Constants.LT_MAX];
		double[] tmpInt2 = new double[Constants.LT_MAX];
                double[] sumInt2 = new double[Constants.LT_MAX];
		double[] tmpSep = new double[Constants.LT_MAX];
		double[] sumSep = new double[Constants.LT_MAX];

		
		double[] sumInt1Sq = new double[Constants.LT_MAX];
	        double[] sumInt2Sq = new double[Constants.LT_MAX];
		double[] sumSepSq = new double[Constants.LT_MAX];
		
		
		
		int[] countSum = new int[Constants.LT_MAX];
		
		TrackDate forSt;
		
		File ferr=null;
		
		int ntr=0, icur=0;
		
		BufferedWriter writer;
		
// Read region file if present

                Region reg = RegionOp.ReadReg();

                if(args.length != 6) {
                  System.err.println("Usage: Det_Mean [path to dirs] [longitude id] [latitude id] [intensity id] [Match dir] [Mean Out Typ. 0 or 1]");
                  System.exit(1);
                }

               
		int ilong = Integer.parseInt(args[1]);
		int ilat = Integer.parseInt(args[2]);
		int ifd = Integer.parseInt(args[3]);
		int ipt = Integer.parseInt(args[5]);
		
		
		if(ipt < 0 || ipt > 1) {
		   System.err.println("Mean output type specifier incorrect, must be 0 or 1");
		   System.exit(1);
		}

                String matchd = args[4];

                File[] forDates = (new File(args[0])).listFiles();

                System.out.println(forDates.length);
		Arrays.sort(forDates);
		
                FilenameFilter fmatchtyp = new FilenameFilter() {
                   public boolean accept(File d, String name) {
                      boolean ret;
                      if(name != null && (name.startsWith("trmatch") )) {
                        ret = true;
                      }
                      else ret = false;
                      return ret;
                   }
                };
		
		for(int i=0; i<Constants.LT_MAX; i++){
								
                        tmpSep[i]=-100;
                        tmpInt1[i]=-100;
                        tmpInt2[i]=-100;
			sumInt1[i]=0;
			sumInt2[i]=0;
			sumSep[i]=0;         
			
			sumInt1Sq[i]=0;
			sumSepSq[i]=0;
			sumInt2Sq[i]=0;
			
			countSum[i]=0;
		
		}
		
		for(int i=0; i<forDates.length; i++){
		
			System.out.println(i + " of " + forDates.length);
		        System.out.println(forDates[i].toString());	
			
			date = forDates[i].toString();
                        date = date.substring(date.length() - 10);
				
			System.out.println(date);
				
			forSt = new TrackDate(Integer.parseInt(date.substring(0,4)), Integer.parseInt(date.substring(4,6)), Integer.parseInt(date.substring(6,8)), Integer.parseInt(date.substring(8,10)));

                        File curdir = new File(forDates[i].toString() + "/" + matchd + "/DET/");

                        if (curdir.isDirectory()) {

                            File[] trackFile = curdir.listFiles(fmatchtyp);

                            if(trackFile != null){
                               //System.out.println(trackFile.length);

                               for(int j=0; j<trackFile.length; j++){
                                  for(int l=0; l<65; l++){
                                     tmpSep[l]=-100;
                                     tmpInt1[l]=-100;
                                     tmpInt2[l]=-100;
                                  }
				
		         	  Track analysis = TrackOperations.File2Track(trackFile[j], 0, ilong, ilat, ifd);//analysis track
				  Track forecast = TrackOperations.File2Track(trackFile[j], 1, ilong, ilat, ifd);//forecast track
				  
// apply region masking to verification track

                                  if(reg.getIreg() == 1) analysis = TrackOperations.TrackMaskRegion(analysis, reg);
				  

			  	  tmpInt2 = intensity(analysis,forecast,forSt, 1);
                                  tmpInt1 = intensity(analysis,forecast,forSt, 0);
				  tmpSep = TrackOperations.separationDistance(analysis,forecast,forSt);
													

				  for(int k=0; k<Constants.LT_MAX; k++){
							
				     if(tmpSep[k] != -100 && tmpInt1[k] != -100 && tmpInt2[k] != -100){
										
				        sumInt1[k]+=tmpInt1[k];
 					sumInt2[k]+=tmpInt2[k];
					sumSep[k]+=tmpSep[k];
					
				        sumInt1Sq[k]+=tmpInt1[k] * tmpInt1[k];
 					sumInt2Sq[k]+=tmpInt2[k] * tmpInt2[k];
					sumSepSq[k]+=tmpSep[k] * tmpSep[k];				
								
					countSum[k]++;
							
				      }

				  }//end of k loop
				  
			       }
			       
			    }
							

			}//end of j for loop
		
		}//end of i for loop


// output deterministic stats

		if(ipt == 0) writer = new BufferedWriter(WriteFile.openfile(null));
		else writer = new BufferedWriter(WriteFile.openfile("det.txt")); 
		
		writer.write("Deterministic Track: Forecast points, Seperation Distance, Intensity Difference, Intensity Bias, Seperation STD, Intensity STD, Intensity Bias STD");
		writer.newLine();
		writer.flush();

                for(int i=0; i<Constants.LT_MAX; i++){
		        sumSep[i]/=countSum[i];
			sumInt1[i]/=countSum[i];
			sumInt2[i]/=countSum[i];
			sumSepSq[i]/=countSum[i];
			sumInt1Sq[i]/=countSum[i];
			sumInt2Sq[i]/=countSum[i];
                        writer.write(Round.round((double)i/4,2) + "\t" + countSum[i] + "\t" 
			                 + Round.round(sumSep[i],4) + "\t" 
					 + Round.round(sumInt1[i],4) + "\t" 
					 + Round.round(sumInt2[i],4) + "\t"
					 + Round.round(Math.sqrt(sumSepSq[i] - sumSep[i] * sumSep[i]),4) + "\t"
					 + Round.round(Math.sqrt(sumInt1Sq[i] - sumInt1[i] * sumInt1[i]),4) + "\t"
					 + Round.round(Math.sqrt(sumInt2Sq[i] - sumInt2[i] * sumInt2[i]),4)       );
		        writer.newLine();
		        writer.flush();
                }
		
		if(ipt == 1) writer.close();
		
	}//end of main method
	
	
        public static double[] intensity(Track analysis, Track forecast, TrackDate forSt, int typ){
	 
           int aPt = 0;
           int fPt = 0;
           TrackDate aDate = new TrackDate(analysis.getPoint(0).getTrackDate());
           TrackDate fDate = new TrackDate(forecast.getPoint(0).getTrackDate());
           TrackDate fstart = new TrackDate(forSt);
		
           boolean flag = true;
           TrackDate max;
           double[] intensity = new double[Constants.LT_MAX];
           for(int i=0; i<Constants.LT_MAX; i++){
	       intensity[i] =-100;
           }

           if(TrackDate.GreaterThan(fDate, aDate)){
	      while(flag){
	          if(TrackDate.GreaterThan(fDate, aDate)){//forecast track starts after analysis track
	    	     aPt ++;
		     aDate.AddTimeFrames(1);
	          }
	          else{
		     flag = false;
	          }
	       }
	    }
		
	    if(TrackDate.GreaterThan(aDate, fDate)){
	      flag = true;	
	      while(flag){
	         if(TrackDate.GreaterThan(aDate, fDate)){//forecast track starts before analysis track
	           fPt ++;
	           fDate.AddTimeFrames(1);
	         }
	         else{
	           flag = false;
	         }
	       }

	    }
		
	    TrackDate aFinish = analysis.getPoint(analysis.getPointNo()-1).getTrackDate();
  	    TrackDate fFinish = forecast.getPoint(forecast.getPointNo()-1).getTrackDate();
		
	    if(TrackDate.GreaterThan(aFinish, fFinish)){
	       max = new TrackDate(fFinish);
	    }
	    else{
	       max =  new TrackDate(aFinish);
	    }
	    int leadTime = 0;
	    flag = true;
	    while(flag){
	       if(TrackDate.Equal(fstart,aDate)){
	          flag = false;
				
	       }
	       else{
	          leadTime++;
	          fstart.AddTimeFrames(1);	
	       }
	    }
		
	    flag = true;
	    int k=0;
		
	    if(!TrackDate.Equal(aDate, fDate)){
	       System.out.println("Error analysis date and forecast date differ");
               System.exit(1);
	    }
			
	    while(flag){
	       if(TrackDate.GreaterThan(aDate,max)){
	          flag = false;
	       }
	       else{
			
	     // System.out.println(forecast.getPoint(fPt) + "\t" + analysis.getPoint(aPt));
	          if(TrackDate.Equal(forecast.getPoint(fPt).getTrackDate(), new TrackDate(0,0,0,0)) ||
	                             (forecast.getPoint(fPt).getLon() > Constants.ADD_CHECK || forecast.getPoint(fPt).getIntensity() > Constants.ADD_CHECK) ||
		                     (analysis.getPoint(aPt).getLon() > Constants.ADD_CHECK || analysis.getPoint(aPt).getIntensity() > Constants.ADD_CHECK)    ){
		     intensity[leadTime] = -100;		 	
	          }
	          else{
		     if(typ == 0){
		        intensity[leadTime] = analysis.getPoint(aPt).getIntensity();
		     }
		     else {
		        intensity[leadTime] = forecast.getPoint(fPt).getIntensity();
		     }
	          }
	          fPt++;
	          aPt++;
	          k++;
	          leadTime++;
	          aDate.AddTimeFrames(1);
	       }

		
	    }
	
	    return intensity;
	
        }
}//end of class
