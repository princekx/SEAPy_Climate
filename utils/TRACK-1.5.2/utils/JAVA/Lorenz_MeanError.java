package tigge;
import java.util.*;
import java.io.*;
import EPSAnalysis.*;

public class Lorenz_MeanError{
	public static void main(String[] args)throws IOException{
		
		
		String date = null;
		String filename = null;
		double[] tmpIntDiff = new double[65];
		double[] sumIntDiff = new double[65];
		double[] tmpSepDiff = new double[65];
		double[] sumSepDiff = new double[65];
		
		double[] sumIntDiff2 = new double[65];
		double[] sumSepDiff2 = new double[65];			
		
		int[] countDiff = new int[65];
		int ntr;
		int trackID;
		
		TrackDate forSt;
		
		if(args.length != 4) {
                  System.err.println("Usage: Lorenz_MeanError [path to dirs] [longitude id] [latitude id] [intensity id]");
                  System.exit(1);
                }
		
                int ilong = Integer.parseInt(args[1]);
		int ilat = Integer.parseInt(args[2]);
		int ifd = Integer.parseInt(args[3]);
		
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
		
		for(int i=0; i<65; i++){
								
                        tmpSepDiff[i]=-100;
                        tmpIntDiff[i]=-100;
			sumIntDiff[i]=0;
			sumSepDiff[i]=0;
			
			sumIntDiff2[i]=0;
			sumSepDiff2[i]=0;
						
			countDiff[i]=0;
		
		}

                File curdir = new File(args[0]);

                if (curdir.isDirectory()) {

                    File[] trackFile = curdir.listFiles(fmatchtyp);
                    Arrays.sort(trackFile);

                    if(trackFile != null){
                       System.out.println(trackFile.length);

                      for(int j=0; j<trackFile.length; j++){
                          for(int l=0; l<65; l++){
                              tmpSepDiff[l]=-100;
                              tmpIntDiff[l]=-100;
                          }
			  
			  System.out.println(trackFile[j].toString());
			  
			  date = trackFile[j].toString().substring(trackFile[j].toString().indexOf("trmatch")).substring(12,22);
			  
			  forSt = new TrackDate(Integer.parseInt(date.substring(0,4)), Integer.parseInt(date.substring(4,6)), Integer.parseInt(date.substring(6,8)), Integer.parseInt(date.substring(8,10)));

                          ntr = TrackOperations.getNoOfTracks(trackFile[j]);
			  
			  for(int k=0; k < ntr; k++){
			      trackID = TrackOperations.getTrackID(trackFile[j], k);
			      
			      if(trackID == 1){
			      
			         Track analysis = TrackOperations.File2Track(trackFile[j], k, ilong, ilat, ifd);//analysis track
			      
			         for(int m=0; m < ntr; m++){
				 
				     trackID = TrackOperations.getTrackID(trackFile[j], m);
				     if(trackID == 2){
				     
				        Track forecast = TrackOperations.File2Track(trackFile[j], m, ilong, ilat, ifd);//forecast track
                                        tmpIntDiff = TrackOperations.intensityDifference(analysis,forecast,forSt);
	                                tmpSepDiff = TrackOperations.separationDistance(analysis,forecast,forSt);
					
				        for(int n=0; n<65; n++){
							
		                           if(tmpSepDiff[n] != -100 && tmpIntDiff[n] != -100){
										
			                      sumIntDiff[n]+=tmpIntDiff[n];
			                      sumSepDiff[n]+=tmpSepDiff[n];
					      
			                      sumIntDiff2[n]+=tmpIntDiff[n] * tmpIntDiff[n];
			                      sumSepDiff2[n]+=tmpSepDiff[n] * tmpSepDiff[n];					      
								
			                      countDiff[n]++;
							
		                           }

		                        }//end of n loop					
							     
				     }	 
				 
				 }
			      
			      }
			  
			  }										
				  
		        }
			       
	              }
							

	          }
		

// output deterministic stats
                System.out.println("Deterministic Track: Forecast points, Seperation Distance, Intensity Difference, Seperation STD, Intensity STD");
                for(int i=0; i<65; i++){
		        sumSepDiff[i]/=countDiff[i];
			sumIntDiff[i]/=countDiff[i];
		        sumSepDiff2[i]/=countDiff[i];
			sumIntDiff2[i]/=countDiff[i];			
                        System.out.println(Round.round((double)i/4,2) + "\t" + countDiff[i] + "\t" 
			                 + Round.round(sumSepDiff[i],4) + "\t" 
					 + Round.round(sumIntDiff[i],4) + "\t"
					 + Round.round(Math.sqrt(sumSepDiff2[i] - sumSepDiff[i] * sumSepDiff[i]),4) + "\t"
					 + Round.round(Math.sqrt(sumIntDiff2[i] - sumIntDiff[i] * sumIntDiff[i]),4)	    );
                }
		
	}//end of main method
}//end of class
