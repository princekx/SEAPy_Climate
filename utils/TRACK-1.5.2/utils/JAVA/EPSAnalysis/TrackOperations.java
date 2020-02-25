/*
 * TrackOperations.java
 *
 * Created on 09 May 2005, 16:08
 */

package EPSAnalysis;

/**
 *
 * @author  lsrf
 */

import java.io.*;
import java.lang.Math.*;
import EPSAnalysis.Constants;

public class TrackOperations {
    
    /** Creates a new instance of TrackOperations */
    public TrackOperations() {
    }

    public static int nthOccurrence(String str, char c, int n) {
       int pos = str.indexOf(c, 0);
       while (n-- > 0 && pos != -1)
           pos = str.indexOf(c, pos+1);
       return pos;
    }
    
    public static Track File2Track(File file, int trackNo, int ilg, int ilt, int ityp)throws IOException, FileNotFoundException{
       
       BufferedReader buffreader = new BufferedReader(new FileReader(file));
       String tmp =  buffreader.readLine();
       boolean flag;
       boolean flag2;
       int ifr = 0;
       int ilng = 0;
       int ilat = 0;
       int pointNo = 0;

       double intensity = 0.0;
       double lon = 0.0;
       double lat = 0.0;

       //System.out.println(file);

       for(int i=0; i<=trackNo; i++){
           flag = true;
           while(flag){
               if(tmp.indexOf("POINT_NUM") == -1){
                   
                   tmp =  buffreader.readLine();
                   //System.out.println(tmp);
               }
               else{
                   flag = false;
               }
           }
			 
           if(i == trackNo){
               
             flag2=true;
	     while(flag2){
		 tmp = tmp.substring(1);
		 if(tmp.substring(0, 1).equals(" ")){
		     flag2 = false;
		 }
	     }
	     flag2=true;
	     while(flag2){
		 tmp = tmp.substring(1);
		 if(!tmp.substring(0, 1).equals(" ")){
	             flag2 = false;
		 }
	     }
	     //System.out.println(tmp);
	     pointNo = Integer.parseInt(tmp);
    
           }
           tmp =  buffreader.readLine();

       }//end of for loop
       
       TrackPoint[] trackpoint = new TrackPoint[pointNo];
		 
       for(int i=0; i<pointNo; i++){
           TrackDate date = new TrackDate(tmp.substring(0,10));
	   if(!TrackDate.Equal(date, new TrackDate(0,0,0,0))){
	   
	        if(ilg == 0){
           	   lon = Double.parseDouble(tmp.substring(11, tmp.indexOf(".")+7));
		}
		else {
		   ilng = nthOccurrence(tmp, '&', ilg-1);
		   lon = Double.parseDouble(tmp.substring(ilng+1, tmp.indexOf("e", ilng+1)+4));
		}
		if(ilt == 0){
           	   lat = Double.parseDouble(tmp.substring(tmp.indexOf(".")+8, tmp.indexOf(".", tmp.indexOf(".")+1)+7));
		}
		else {
		   ilat = nthOccurrence(tmp, '&', ilt-1);
		   lat = Double.parseDouble(tmp.substring(ilat+1, tmp.indexOf("e", ilat+1)+4));
		}
                if(ityp == 0) {
           	   intensity = Double.parseDouble(tmp.substring(tmp.indexOf(".", tmp.indexOf(".")+1)+8, tmp.indexOf("e")+4));
                }
                else {
		   ifr = nthOccurrence(tmp, '&', ityp-1);
                   intensity = Double.parseDouble(tmp.substring(ifr+1, tmp.indexOf("e", ifr+1)+4)); 
                }
		
		//System.out.println(lon + "\t" + lat + "\t" + intensity);
		//System.exit(1);
		
                trackpoint[i] = new TrackPoint(lon, lat, date, intensity);
	   }
	   else{
		trackpoint[i] = new TrackPoint(0,0,date,0);
	   }
	   
           if(i!=(pointNo-1)){ 
                tmp =  buffreader.readLine();
           }
       }
       
       Track track = new Track(trackpoint);
		 
		 buffreader.close();
       return track;
       
    }
	 
    public static int getTrackID(File file, int trackNo)throws IOException, FileNotFoundException{
       
       BufferedReader buffreader = new BufferedReader(new FileReader(file));
       String tmp =  buffreader.readLine();
       boolean flag;
		 boolean flag2;
       int pointNo = 0;
       for(int i=0; i<=trackNo; i++){
           flag = true;
           while(flag){
			 
               if(tmp.indexOf("TRACK_ID") == -1){
                   
                   tmp =  buffreader.readLine();
               }
               else{
                   flag = false;
               }
           }
			
          if(i!=trackNo){
                tmp =  buffreader.readLine();
          }

       }//end of for loop
		 
       tmp = tmp.substring(8);	
       flag = true;
       while(flag){
	  if(tmp.indexOf(" ")==0){
	     tmp = tmp.substring(1);
	  }
	  else{
             flag = false;	
				
	  }
       }
       int trackID = Integer.parseInt(tmp.substring(0, tmp.indexOf(" ")));
       buffreader.close();
       return trackID;
			
    } 
    public static int getNoOfTracks(File file)throws IOException, FileNotFoundException{
	 	BufferedReader buffreader = new BufferedReader(new FileReader(file));
		String tmp = buffreader.readLine();
		boolean flag;
		
	        flag = true;
                while(flag){
                   if(tmp.indexOf("TRACK_NUM") == -1){
                   
                       tmp =  buffreader.readLine();
                   }
                   else{
                       flag = false;
                   }
                }
//		tmp = buffreader.readLine();
//		tmp = buffreader.readLine();
		tmp = tmp.substring(9);
		flag=true;
		while(flag){
			if(tmp.indexOf(" ")==0){
				tmp = tmp.substring(1);
			}
			else{
				flag = false;
				if(tmp.indexOf(" ") != -1){
					tmp = tmp.substring(0, tmp.indexOf(" "));
				}
			}
		}
	 	int noOfTracks = Integer.parseInt(tmp);
		buffreader.close();
		return noOfTracks;
	 
	 }
	 
    public static double[] intensityDifference(Track analysis, Track forecast, TrackDate forSt){
	 
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
		 intensity[leadTime] = Math.abs(forecast.getPoint(fPt).getIntensity() - analysis.getPoint(aPt).getIntensity());
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
	 
	 
	 
    public static double[] separationDistance(Track analysis, Track forecast, TrackDate forSt){
	 	
	int aPt = 0;
	int fPt = 0;
	TrackDate aDate = new TrackDate(analysis.getPoint(0).getTrackDate());
	TrackDate fDate = new TrackDate(forecast.getPoint(0).getTrackDate());
	TrackDate fstart = new TrackDate(forSt);
        TrackDate zDate = new TrackDate(0,0,0,0);
		
	boolean flag = true;
	TrackDate max;
	double[] sep = new double[Constants.LT_MAX];
	
	for(int i=0; i<Constants.LT_MAX; i++){
	    sep[i] =-100;
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
		
	int z=0;
	while(flag){
				
	   if(TrackDate.Equal(fstart,aDate)){
	      flag = false;
				
	   }
	   else{
	      leadTime++;
	      fstart.AddTimeFrames(1);	
	   }
	   z++;
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
		if(TrackDate.Equal(forecast.getPoint(fPt).getTrackDate(), zDate)                  ||
                   TrackDate.Equal(analysis.getPoint(aPt).getTrackDate(), zDate)                  ||
			           forecast.getPoint(fPt).getLon() > Constants.ADD_CHECK          ||
			           analysis.getPoint(aPt).getLon() > Constants.ADD_CHECK             ){
		   sep[leadTime] = -100;
				
		}
		else{
		   sep[leadTime] = geodetic(forecast.getPoint(fPt), analysis.getPoint(aPt));
		}
		fPt++;
		aPt++;
		k++;
		leadTime++;
		aDate.AddTimeFrames(1);
	     }

		
	}
	
	return sep;
	
    }
	 
	 
	 
	 public static double MeanSeparation(Track track1, Track track2, TrackDate forSt){//calculates mean separation distance of two tracks
	 	
		double[] tmpSep = new double[Constants.LT_MAX];
		tmpSep = separationDistance(track1, track2, forSt);
		
		double meanSep = 0.0;
		int count = 0;
		
		for(int i=0; i<Constants.LT_MAX; i++){
	 		if(tmpSep[i]!=-100){
				meanSep += tmpSep[i];
				count ++;
				
			}
	 
	 	}
		
		meanSep = meanSep/count;
		return meanSep;
	 
	 }
	 
	 public static double MeanSeparation(Track track1, Track track2, TrackDate forSt, int min, int max){//calculates mean separation distance of two tracks between min and max
	 	
		double[] tmpSep = new double[Constants.LT_MAX];
		tmpSep = separationDistance(track1, track2, forSt);
		
		double meanSep = 0.0;
		int count = 0;
		
		for(int i=min; i<max; i++){
	 		if(tmpSep[i]!=-100){
				meanSep += tmpSep[i];
				count ++;
				
			}
	 
	 	}
		
		meanSep = meanSep/count;
		return meanSep;
	 
	 }
	 
	 
	 
	 public static double MeanIntensity(Track track1, Track track2, TrackDate forSt){//calculates mean absolute intensty difference of two tracks
	 	
		double[] tmpInt = new double[Constants.LT_MAX];
		tmpInt = intensityDifference(track1, track2, forSt);
		
		double meanInt = 0.0;
		int count = 0;
		
		for(int i=0; i<Constants.LT_MAX; i++){
	 		if(tmpInt[i]!=-100){
				meanInt += tmpInt[i];
				count ++;
				
			}
	 
	 	}
		
		meanInt = meanInt/count;
		return meanInt;
	 
	 }
	 
	 
	 public static double MeanIntensity(Track track1, Track track2, TrackDate forSt, int min, int max){//calculates mean absolute intensty difference of two tracks
	 	
		double[] tmpInt = new double[Constants.LT_MAX];
		tmpInt = intensityDifference(track1, track2, forSt);
		
		double meanInt = 0.0;
		int count = 0;
		
		for(int i=min; i<max; i++){
	 		if(tmpInt[i]!=-100){
				meanInt += tmpInt[i];
				count ++;
				
			}
	 
	 	}
		
		meanInt = meanInt/count;
		return meanInt;
	 
	 }
	 
    public static double geodetic(TrackPoint point1, TrackPoint point2){
	 	
	if((point1.getLon() == point2.getLon()) && (point1.getLat() == point2.getLat())){
	   return 0.0;	
	}
		
	double lon1 = point1.getLon()*Math.PI/180;//point 1 longitude in radians
	double lon2 = point2.getLon()*Math.PI/180;//point 2 longitude in radians
		
	double lat1 = point1.getLat()*Math.PI/180;//point 1 latitude in radians
	double lat2 = point2.getLat()*Math.PI/180;//point 1 latitude in radians
		
	double dist = Math.acos(Math.sin(lat1)*Math.sin(lat2)+Math.cos(lat1)*Math.cos(lat2)*Math.cos(lon2-lon1));
		
	return dist*180/Math.PI;
	 	 
    }	 
	 
//actual intensity difference, not absoloute
	 
    public static double[] intensityDifference2(Track analysis, Track forecast, TrackDate forSt){
	 
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
			
		//System.out.println(forecast.getPoint(fPt) + "\t" + analysis.getPoint(aPt));
		if(TrackDate.Equal(forecast.getPoint(fPt).getTrackDate(), new TrackDate(0,0,0,0)) ||
	                           (forecast.getPoint(fPt).getLon() > Constants.ADD_CHECK || forecast.getPoint(fPt).getIntensity() > Constants.ADD_CHECK) ||
		                   (analysis.getPoint(aPt).getLon() > Constants.ADD_CHECK || analysis.getPoint(aPt).getIntensity() > Constants.ADD_CHECK)     ){
	  	   intensity[leadTime] = -100;
				
		}
		else{
		   intensity[leadTime] = forecast.getPoint(fPt).getIntensity() - analysis.getPoint(aPt).getIntensity();
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

	 //calculates percentage intensity difference
	 public static double[] intensityDifferencePercentage(Track analysis, Track forecast, TrackDate forSt){
	 
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
			
				//System.out.println(forecast.getPoint(fPt) + "\t" + analysis.getPoint(aPt));
				if(TrackDate.Equal(forecast.getPoint(fPt).getTrackDate(), new TrackDate(0,0,0,0))){
					intensity[leadTime] = -100;
				
				}
				else{
					intensity[leadTime] = 100*(Math.abs(forecast.getPoint(fPt).getIntensity() - analysis.getPoint(aPt).getIntensity()))/analysis.getPoint(aPt).getIntensity();
					/*if(intensity[leadTime]>100){
						System.out.println("analysis");
						System.out.println("leadTime = " + leadTime);
						
						System.out.println(analysis + "\n");
						System.out.println(forecast + "\n");
					}*/
				
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
	 
	 
	 
	 	//actual intensity difference, not absoloute	 
	 	public static double[] analIntLeadTime(Track analysis, Track forecast, TrackDate forSt){
	 
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
			
				//System.out.println(forecast.getPoint(fPt) + "\t" + analysis.getPoint(aPt));
				if(TrackDate.Equal(forecast.getPoint(fPt).getTrackDate(), new TrackDate(0,0,0,0))){
					intensity[leadTime] = -100;
				
				}
				else{
					intensity[leadTime] = analysis.getPoint(aPt).getIntensity();
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

	 public static boolean region(Track track, double minLat, double maxLat){// method to check if track lies in latitude band, returns true if it does
	 
	 	boolean flag = true;
		int count = 0;
		
		while(flag && (count < track.getPointNo()) ){
			
	 		if((track.getPoint(count).getLat() > minLat) && (track.getPoint(count).getLat() < maxLat)){
				count++;
			
			}
			else{
				flag = false;
			
			}
	 
	 
	 	}//end of while
	 
	 	return flag;
	 }//end of region method

         public static Track TrackMaskRegion(Track track, Region reg) { // method to mask track by region

               TrackPoint trp;

               for(int i=0; i < track.getPointNo(); i++){
                  trp = track.getPoint(i);
                  if(RegionOp.checkRegion(reg, trp) == 0){
                     trp.setLon(Constants.ADD_UNDEF);
                     trp.setLat(Constants.ADD_UNDEF);
                     trp.setIntensity(Constants.ADD_UNDEF);
                  }
               } 

               return track;
         }
	 
	 public static Track TrackMaskTime(Track track1, Track track2){ // method to mask one track by another that has already been region filtered.
	 
               int[] maskTrack2 = new int[track2.getPointNo()];

	       TrackPoint trp1, trp2;
	       
//	       System.out.println("+++++++");
//	       for(int i=0; i < track1.getPointNo(); i++){
//	           trp1 = track1.getPoint(i);
//	           System.out.println(trp1.getTrackDate() + "\t" + trp1.getLon() + "\t" + trp1.getLat());
//	       }
	       
	       for(int i=0; i < track1.getPointNo(); i++){
	           trp1 = track1.getPoint(i);
		   for(int j=0; j < track2.getPointNo(); j++){
		      trp2 = track2.getPoint(j);
		      if(TrackDate.Equal(trp1.getTrackDate(), trp2.getTrackDate())){
		         maskTrack2[j] = 1;
		         if(trp1.getLon() > Constants.ADD_CHECK){
			    trp2.setLon(Constants.ADD_UNDEF);
                            trp2.setLat(Constants.ADD_UNDEF);
                            trp2.setIntensity(Constants.ADD_UNDEF);
			 }
			 
		      }
		        		      
		   }
	       
	       }
	       for(int j=0; j < track2.getPointNo(); j++){
	           trp2 = track2.getPoint(j);
	           if(maskTrack2[j] == 0){
		       trp2.setLon(Constants.ADD_UNDEF);
                       trp2.setLat(Constants.ADD_UNDEF);
                       trp2.setIntensity(Constants.ADD_UNDEF); 
		   }
	       }
	       
//	       System.out.println(".......");
//	       for(int i=0; i < track2.getPointNo(); i++){
//	           trp2 = track2.getPoint(i);
//	           System.out.println(trp2.getTrackDate() + "\t" + trp2.getLon() + "\t" + trp2.getLat());
//	       }
	 
	       return track2;
	 }
    
}
