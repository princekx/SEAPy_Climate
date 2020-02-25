package EPSAnalysis;

public class TrackPoint{

	private double lon;
	private double lat;
	private TrackDate trackdate;
	private double intensity;
	
	TrackPoint(double lon, double lat, TrackDate trackdate, double intensity){
	
	this.lon = lon;
	this.lat = lat;
	this.trackdate = trackdate;
	this.intensity = intensity;
	
	}
	
	TrackPoint(TrackPoint trackpoint){
		
		this.lon = trackpoint.getLon();
		this.lat = trackpoint.getLat();
		this.trackdate = trackpoint.getTrackDate();
		this.intensity = trackpoint.getIntensity();
		
	
	}
	
	
	public void setLon(double lon){
	
		this.lon = lon;
	
	}	
	
	public void setLat(double lat){
	
		this.lat = lat;
	
	}	
	
	public void setTrackDate(TrackDate trackdate){
	
		this.trackdate = trackdate;
	
	}	
	
	public void setIntensity(double intensity){
	
		this.intensity = intensity;
	
	}	

	public double getLon(){
	
		return lon;
	
	}	
	
	public double getLat(){
	
		return lat;
	
	}	
	
	public TrackDate getTrackDate(){
	
		return trackdate;
	
	}	
	
	public double getIntensity(){
	
		return intensity;
	
	}	


	public String toString(){

	return trackdate + " " + lon + " " + lat + " " + intensity;
	


	}


}
