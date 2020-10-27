# VelodromeModel

VelodromeModel creates a velodrome track black-line model that consists of two straights, two circular arc bends and four transition curves between the bends and the straights. The transition curves follow the form of a clothoid (also known as a Euler spiral or Cornu spiral). A clothoid has the unique property that its curvature is a polynomial (usually linear) function of its arc length. This provides controlling of the centripetal acceleration during cornering and has G2 geometric continuity.

The features that define the track are:
   * L_L: The lap length. This is generally a known, fixed value. 
   * Y: The half-width between the two straights.
   * R: The turn radius at the bend apex.
   * n: The power of the change in curvature with length. Usually 1 (linear)
   
The ratio of inputs Y/R has a limited range of feasible solutions dependent on L_L, R, & n that are checked before calculations begin. 
