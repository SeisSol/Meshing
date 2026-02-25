#include "SeismicVelocity.h"

double landers61(int, double, double, double z)
{
  if (z > -0.1) {
    return 1925.9102873142;
  }
  if (z >= -0.3) {
    return 3473.1266454957;
  }
  if (z >= -1.0) {
    return 4224.2887430448;
  }
  if (z >= -3.0) {
    return 5100.0000000000;
  }
  if (z >= -6.0) {
    return 5922.6430205131;
  }
  if (z >= -31.0) {
    return 5933.4162648030;
  } else {
    return 7800.0000000000;
  }
  return -1.0;
}

/*
BedrockVelModel(1,:) = (/  -6d3, 2550d0,18589500000d0,26571000000d0/) => 5000 m/s
BedrockVelModel(2,:) = (/  -8d3, 2850d0,39016500000d0,42379500000d0/) => 6500 m/s
BedrockVelModel(3,:) = (/ -12d3, 3050d0,50027625000d0,53695250000d0/) => 7100 m/s
BedrockVelModel(4,:) = (/-6d3,2720d0,33320000000d0,31280000000d0/)    => 6000 m/s
BedrockVelModel(5,:) = (/-12d3,2860d0,41298400000d0,41984800000d0/)   => 6600 m/s
BedrockVelModel(6,:) = (/-23d3,3050d0,46390500000d0,60969500000d0/)   => 7100 m/s
BedrockVelModel(7,:) = (/ -5d10, 3330d0,65942325000d0,81235350000d0/) => 8000 m/s  */
double sumatra1223_high(int group, double, double, double z)
{
  double cp = -1.0;
  switch (group) {
    case 4: // LVZ
    case 2: // layer 1
      cp = 5000.0;
      break;
    case 3: // layer 2
      cp = 6500.0;
      break;
    case 7: // layer 3
      cp = 7100.0;
      break;
    case 6: // layer 4
      cp = 8000.0;
      break;
    case 1: // continental
    case 5: // box
      if (z > -6000.0) {
        cp = 6000.0;
      } else if (z >= -12000.0) {
        cp = 6600.0;
      } else if (z > -23000.0) {
        cp = 7100.0;
      } else {
        cp = 8000.0;
      }
      break;
    default:
      break;
  }
  return cp;
}

double sumatra1223_low(int group, double, double, double z)
{
  double cp = -1.0;
  switch (group) {
    case 2: // LVZ
    case 7: // layer 1
      cp = 5000.0;
      break;
    case 4: // layer 2
      cp = 6500.0;
      break;
    case 3: // layer 3
      cp = 7100.0;
      break;
    case 6: // layer 4
      cp = 8000.0;
      break;
    case 1: // continental
    case 5: // box
      if (z > -6000.0) {
        cp = 6000.0;
      } else if (z >= -12000.0) {
        cp = 6600.0;
      } else if (z > -23000.0) {
        cp = 7100.0;
      } else {
        cp = 8000.0;
      }
      break;
    default:
      break;
  }
  return cp;
}

double sumatra1224(int group, double, double, double z)
{
  double cp = -1.0;
  switch (group) {
    case 1:
      cp = 6500.0;
      break;
    case 2:
      cp = 5000.0;
      break;
    case 3:
      cp = 7100.0;
      break;
    case 4:
    case 5:
      cp = 8000.0;
      break;
    case 6:
      if (z > -6000.0) {
        cp = 6000.0;
      } else if (z >= -12000.0) {
        cp = 6600.0;
      } else if (z > -23000.0) {
        cp = 7100.0;
      } else {
        cp = 8000.0;
      }
      break;
    default:
      break;
  }
  return cp;
}
