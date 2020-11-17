# VelodromeModel

VelodromeModel creates a velodrome track black-line model that consists of two straights, two circular arc bends and four transition curves between the bends and the straights. The transition curves are based on two different Cesaro equations where the curvature along the arc length is defined. This provides controlling of the centripetal acceleration during cornering and has either G2 or G3 geometric continuity. 

The first option is for a generalised clothoid (also known as a Euler spiral or Cornu spiral). Here the curvature of the transition curve is a polynomial function of its arc length and the exponent power, _n_, is required to be selected. A typical case is for _n_ = 1 which is the standard Euler spiral transition curve. This is a G2 continuous curve. 

The second option is for a half-sine wave curvature profile where the curvature increases from zero to the bend curvature following a sinusoidal path. This is a G3 continuous curve. 

The measurable features that define the track are:
   * L_L: The lap length. This is generally a known, fixed value. 
   * Y: The half-width between the two straights.
   * R: The turn radius at the bend apex.
   
The ratio of inputs Y/R has a limited range of feasible solutions dependent on L_L, R, and the curvature function that are checked before calculations begin. 
