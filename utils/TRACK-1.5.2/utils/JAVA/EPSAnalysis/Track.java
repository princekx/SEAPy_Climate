package EPSAnalysis;

public class Track{

	
	private int noOfPoints = 0;
	
	private TrackPoint[] trackpoints;
	
	Track(TrackPoint[] trackpoints){
		
		noOfPoints = trackpoints.length;
		this.trackpoints = new TrackPoint[noOfPoints];
		
		for(int i=0; i< noOfPoints; i++){
		
			this.trackpoints[i] = trackpoints[i];
		}
	
	}
        
        Track(Track track){
            
            this.noOfPoints = track.getPointNo();
            trackpoints = new TrackPoint[noOfPoints];
            
            for(int i=0; i<noOfPoints; i++){
                this.trackpoints[i] = track.getPoint(i);
                
            }
            
        }
        
        public TrackPoint getPoint(int pointNo){
            if(pointNo >= noOfPoints){
               System.out.println("Point number to large!!"); 
            }
            return trackpoints[pointNo];
        }
        
        public int getPointNo(){
            return noOfPoints;
        }
		  
		  public void setPointNo(int noOfPoints){
		  		this.noOfPoints = noOfPoints;
		  }
        public void setPoint(int pointNo, TrackPoint trackpoint){
            trackpoints[pointNo] = new TrackPoint(trackpoint);
        }
	
	public String toString(){
            
            String tmp = "";
            for(int i=0; i<noOfPoints; i++){
                tmp += (trackpoints[i] + "\n");
               
            }
            return tmp;
            
  }
  
  public double getMaxIntensity(){
  		
		double maxIntensity =0;
		for(int i=0; i<trackpoints.length; i++){
  
  			if(getPoint(i).getIntensity() > maxIntensity){
				maxIntensity = getPoint(i).getIntensity();
			}
  
  
  		}
		
		return maxIntensity;
  
  }


}
