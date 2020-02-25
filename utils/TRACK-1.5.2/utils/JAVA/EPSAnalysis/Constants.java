/*
 * Constants.java
 *
 */

package EPSAnalysis;

/**
 *
 * @author  kih
 */

public final class Constants {

/** Array dimensions */

  public static final int LT_MAX = 65;

/** Constants for checking missing values */

  public static final Double ADD_UNDEF = 1.0e+25;
  public static final Double ADD_CHECK = 1.0e+20;

  private Constants(){
    //this prevents even the native class from 
    //calling this ctor as well :
    throw new AssertionError();
  }

}
