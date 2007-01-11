#include "Isosurface.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/gle.h>


// EdgeData Isosurface::CubeEdgeData[256] = {      /// 76543210
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00000000
//   { 1,    1, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00000001
//   { 1,    1, 6,11, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00000010
//   { 2,    5, 6,11, 5, 9,11, 0, 0, 0, 0, 0, 0 }, /// 00000011
//   { 1,    5,10,11, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00000100
//   { 2,    1, 9,10, 1,11, 9, 0, 0, 0, 0, 0, 0 }, /// 00000101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00000110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00000111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00001000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00001001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00001010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00001011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00001100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00001101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00001110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00001111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00010000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00010001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00010010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00010011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00010100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00010101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00010110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00010111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00011000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00011001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00011010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00011011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00011100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00011101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00011110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00011111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00100000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00100001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00100010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00100011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00100100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00100101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00100110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00100111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00101000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00101001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00101010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00101011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00101100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00101101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00101110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00101111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00110000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00110001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00110010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00110011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00110100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00110101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00110110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00110111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00111000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00111001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00111010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00111011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00111100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00111101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00111110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 00111111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01000000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01000001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01000010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01000011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01000100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01000101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01000110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01000111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01001000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01001001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01001010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01001011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01001100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01001101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01001110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01001111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01010000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01010001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01010010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01010011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01010100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01010101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01010110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01010111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01011000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01011001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01011010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01011011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01011100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01011101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01011110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01011111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01100000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01100001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01100010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01100011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01100100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01100101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01100110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01100111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01101000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01101001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01101010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01101011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01101100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01101101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01101110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01101111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01110000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01110001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01110010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01110011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01110100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01110101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01110110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01110111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01111000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01111001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01111010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01111011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01111100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01111101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01111110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 01111111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10000000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10000001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10000010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10000011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10000100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10000101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10000110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10000111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10001000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10001001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10001010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10001011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10001100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10001101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10001110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10001111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10010000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10010001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10010010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10010011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10010100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10010101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10010110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10010111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10011000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10011001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10011010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10011011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10011100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10011101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10011110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10011111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10100000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10100001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10100010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10100011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10100100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10100101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10100110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10100111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10101000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10101001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10101010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10101011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10101100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10101101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10101110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10101111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10110000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10110001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10110010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10110011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10110100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10110101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10110110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10110111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10111000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10111001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10111010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10111011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10111100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10111101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10111110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 10111111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11000000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11000001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11000010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11000011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11000100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11000101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11000110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11000111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11001000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11001001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11001010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11001011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11001100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11001101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11001110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11001111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11010000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11010001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11010010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11010011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11010100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11010101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11010110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11010111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11011000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11011001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11011010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11011011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11011100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11011101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11011110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11011111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11100000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11100001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11100010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11100011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11100100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11100101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11100110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11100111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11101000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11101001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11101010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11101011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11101100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11101101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11101110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11101111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11110000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11110001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11110010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11110011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11110100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11110101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11110110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11110111
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11111000
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11111001
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11111010
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11111011
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11111100
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11111101
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// 11111110
//   { 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  /// 11111111
// };


int
Isosurface::CubicFormula (double a, double b, double c, double d,
			  double &x1, double &x2, double &x3)
{
  double A = b/a;
  double B = c/a;
  double C = d/a;
  double Q = (A*A - 3.0*B)/9.0;
  double R = (2.0*A*A*A - 9.0*A*B + 27.0*C)/54.0;
  const double third = 1.0/3.0;
  const double thirdA = third * A;
  //cerr << "Q = " << Q << " R = " << R << "\n";
  if ((R*R) < (Q*Q*Q)) {
    double theta = acos(R/sqrt(Q*Q*Q));
    double twosqrtQ = 2.0*sqrt(Q);
    x1 = -twosqrtQ*cos(third*theta) - thirdA;
    x2 = -twosqrtQ*cos(third*(theta + 2.0*M_PI)) - thirdA;
    x3 = -twosqrtQ*cos(third*(theta - 2.0*M_PI)) - thirdA;
    return 3;
  }
  else {
    double sign = (R>0.0) ? 1.0 : -1.0;
    double alpha = -sign*cbrt(fabs(R) + sqrt((R*R)-(Q*Q*Q)));
    double beta = (alpha == 0.0) ? 0.0 : (Q/alpha);
    x1 = (alpha + beta) - thirdA;
    double check = a*(x1*x1*x1)+b*(x1*x1)+c*x1+d;
    double maxcoef = max(fabs(a),max(fabs(b),max(fabs(c),fabs(d))));
    assert (fabs(check) < (1e-12*maxcoef));
    return 1;
  }
}



Vec3
Isosurface::FindEdge (int ix, int iy, int iz, int edgeNum)
{
  int ix1, ix2, iy1, iy2, iz1, iz2;
  int dim;
  dim = EdgeTable[edgeNum][0];
  ix1 = ix+EdgeTable[edgeNum][1];   ix2 = ix+EdgeTable[edgeNum][2];
  iy1 = iy+EdgeTable[edgeNum][3];   iy2 = iy+EdgeTable[edgeNum][4];
  iz1 = iz+EdgeTable[edgeNum][5];   iz2 = iz+EdgeTable[edgeNum][6];

  Vec3 r1((*Xgrid)(ix1),(*Ygrid)(iy1),(*Zgrid)(iz1));
  Vec3 r2((*Xgrid)(ix2),(*Ygrid)(iy2),(*Zgrid)(iz2));
  Vec3 delta  = r2-r1;
  
  double v1  = F(ix1,iy1,iz1)[0] - Isoval; 
  double v2  = F(ix2,iy2,iz2)[0] - Isoval;

  assert ((v1*v2) < 0.0);
  double u;
  if (UseCubicInterp) {
    double dvdu1 = F(ix1,iy1,iz1)[dim+1] * delta[dim];
    double dvdu2 = F(ix2,iy2,iz2)[dim+1] * delta[dim];
    double a, b, c, d, u1, u2, u3;
    a = dvdu1 + dvdu2 + 2.0*(v1-v2);
    b = 3.0*(v2-v1) - 2.0*dvdu1 - dvdu2;
    c = dvdu1;
    d = v1;
    int numSols = CubicFormula (a, b, c, d, u1, u2, u3);
    if (numSols == 1)
      u = u1;
    else if ((u1>=0.0) && (u1<=1.0))
      u = u1;
    else if ((u2>=0.0) && (u2<=1.0))
      u = u2;
    else 
      u = u3;
    assert ((u>=0.0) && (u<=1.0));
  }
  else   // Linear interpolation
    u = fabs(v1/(v2-v1));
  
  return r1 + u*delta;
}

void
Isosurface::Set()
{
  if (!UpToDate)
    Update();
  Start();

  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30.0);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  if (UseNormals)
    glEnable (GL_NORMALIZE);

  glDepthMask(GL_FALSE);

  int numTriangles = 0;

  for (int i=0; i < Isovals.size(); i++) {
    Isoval = Isovals[i];
    Color  = Colors[i];
    glColor4d(Color[0], Color[1], Color[2], Color[3]);
    float fcolor[4] = { Color[0], Color[1], Color[2], Color[3] };
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
    float spec[4] = { 1.0, 1.0, 1.0, Color[3]};
    glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, spec);

    glBegin(GL_TRIANGLES);
    for (int ix=0; ix<(Nx-1); ix++) {
      for (int iy=0; iy<(Ny-1); iy++) {
	for (int iz=0; iz<(Nz-1); iz++) { 
	  /// Check corners
	  int index = 0;
	  index |= (F(ix  ,  iy+1, iz  )[0] > Isoval);
	  index |= ((F(ix+1, iy+1, iz  )[0] > Isoval) << 1);
	  index |= ((F(ix+1, iy  , iz  )[0] > Isoval) << 2);
	  index |= ((F(ix  , iy  , iz  )[0] > Isoval) << 3);
	  index |= ((F(ix  , iy+1, iz+1)[0] > Isoval) << 4);
	  index |= ((F(ix+1, iy+1, iz+1)[0] > Isoval) << 5);
	  index |= ((F(ix+1, iy  , iz+1)[0] > Isoval) << 6);
	  index |= ((F(ix  , iy  , iz+1)[0] > Isoval) << 7);
	  int ei=0;
	  int edge;
	  while ((edge=EdgeData[index][ei]) != -1) {
	    numTriangles++;
	    Vec3 reduced = FindEdge (ix, iy, iz, edge);
// 	    Vec3 vertex = (reduced[0]*LatticeVecs[0] +
// 			   reduced[1]*LatticeVecs[1] +
// 			   reduced[2]*LatticeVecs[2]);
	    Vec3 vertex = reduced * Lattice;
	    if (UseNormals) {
	      Vec3 nred = -1.0*Grad(reduced[0], reduced[1], reduced[2]);
// 	      Vec3 normal = (nred[0]*LatticeVecs[0]+
// 			     nred[1]*LatticeVecs[1]+
// 			     nred[2]*LatticeVecs[2]);
	      Vec3 normal = nred*Lattice;
	      glNormal3dv(&(normal[0]));
	    }
	    glVertex3dv(&(vertex[0]));
	    ei++;
	  }
	}
      }
    }
    glEnd();
  }
  glDepthMask(GL_TRUE);
  End();
}


int
Isosurface::NumTriangles(int i)
{
  int num = 0;
  Isoval = Isovals[i];
  for (int ix=0; ix<(Nx-1); ix++) {
    for (int iy=0; iy<(Ny-1); iy++) {
      for (int iz=0; iz<(Nz-1); iz++) { 
	/// Check corners
	int index = 0;
	index |= (F(ix  ,  iy+1, iz  )[0] > Isoval);
	index |= ((F(ix+1, iy+1, iz  )[0] > Isoval) << 1);
	index |= ((F(ix+1, iy  , iz  )[0] > Isoval) << 2);
	index |= ((F(ix  , iy  , iz  )[0] > Isoval) << 3);
	index |= ((F(ix  , iy+1, iz+1)[0] > Isoval) << 4);
	index |= ((F(ix+1, iy+1, iz+1)[0] > Isoval) << 5);
	index |= ((F(ix+1, iy  , iz+1)[0] > Isoval) << 6);
	index |= ((F(ix  , iy  , iz+1)[0] > Isoval) << 7);
	int ei=0;
	while (EdgeData[index][ei] != -1) {
	  num++;
	  ei++;
	}
      }
    }
  }
  return num/3;
}

void
Isosurface::DrawPOV (FILE *fout, string rotString)
{
  for (int i=0; i < Isovals.size(); i++) {
    if (NumTriangles(i) != 0) {
      Isoval = Isovals[i];
      Color  = Colors[i];
      fprintf (fout, "mesh {\n");
      for (int ix=0; ix<(Nx-1); ix++) {
	for (int iy=0; iy<(Ny-1); iy++) {
	  for (int iz=0; iz<(Nz-1); iz++) { 
	    /// Check corners
	    int index = 0;
	    index |= (F(ix  ,  iy+1, iz  )[0] > Isoval);
	    index |= ((F(ix+1, iy+1, iz  )[0] > Isoval) << 1);
	    index |= ((F(ix+1, iy  , iz  )[0] > Isoval) << 2);
	    index |= ((F(ix  , iy  , iz  )[0] > Isoval) << 3);
	    index |= ((F(ix  , iy+1, iz+1)[0] > Isoval) << 4);
	    index |= ((F(ix+1, iy+1, iz+1)[0] > Isoval) << 5);
	    index |= ((F(ix+1, iy  , iz+1)[0] > Isoval) << 6);
	    index |= ((F(ix  , iy  , iz+1)[0] > Isoval) << 7);
	    int ei=0;
	    int edge;
	    int triCounter = 0;
	    while ((edge=EdgeData[index][ei]) != -1) {
	      if (triCounter == 0)
		fprintf (fout, "  smooth_triangle {\n");
	      
	      Vec3 reduced = FindEdge (ix, iy, iz, edge);
	      Vec3 vertex = (reduced[0]*LatticeVecs[0] +
			     reduced[1]*LatticeVecs[1] +
			     reduced[2]*LatticeVecs[2]);
	      Vec3 nred = -1.0*Grad(reduced[0], reduced[1], reduced[2]);
	      Vec3 normal = (nred[0]*LatticeVecs[0]+
			     nred[1]*LatticeVecs[1]+
			     nred[2]*LatticeVecs[2]);
	      normal = 1.0/sqrt(dot(normal,normal)) * normal;
	      fprintf (fout, "    <%14.10f, %14.10f, %14.10f>, ",
		       vertex[0], vertex[1], vertex[2]);
	      fprintf (fout, " <%14.10f, %14.10f, %14.10f>",
		       normal[0], normal[1], normal[2]);
	      
	      triCounter++;
	      if (triCounter == 3) {
		fprintf (fout, "\n  }\n");
		triCounter = 0;
	      }
	      else
		fprintf (fout, ",\n");
	      ei++;
	    }
	  }
	}
      }
      fprintf (fout, "  pigment { color rgbt <%1.5f %1.5f %1.5f %1.5f> }\n", 
	       0.6*Color[0], 0.6*Color[1], 0.6*Color[2], 1.0-0.7*Color[3]);
      fprintf (fout, "  finish { \n");
      fprintf (fout, "    specular 0.25roughness 0.025\n");
      fprintf (fout, "    ambient  0.2\n");
      fprintf (fout, "    diffuse  0.8\n");
      fprintf (fout, "  }\n");
      fprintf (fout, "%s", rotString.c_str());
      fprintf (fout, "}\n\n");
    }
  }
}


void
Isosurface::SetColor (TinyVector<double,3> color)
{
  Colors.resize(1);
  Colors[0][0] = color[0];
  Colors[0][1] = color[1];
  Colors[0][2] = color[2];
  Colors[0][3] = Alpha;
}

void
Isosurface::SetAlpha(double alpha)
{
  Alpha = alpha;
  for (int i=0; i<Colors.size(); i++)
    Colors[i][3] = alpha;
}


void
Isosurface::SetColor (vector<TinyVector<double,3> > &colors)
{
  Colors.resize(colors.size());
  for (int i=0; i<colors.size(); i++) {
    Colors[i][0] = colors[i][0];
    Colors[i][1] = colors[i][1];
    Colors[i][2] = colors[i][2];
    Colors[i][3] = Alpha;
  }
}


int Isosurface::EdgeTable[12][7] = {
// dim  x1 x2 y1 y2 z1 z2
  { 0,   0, 1, 1, 1, 0, 0},
  { 1,   1, 1, 0, 1, 0, 0},
  { 0,   0, 1, 0, 0, 0, 0},
  { 1,   0, 0, 0, 1, 0, 0},
  { 0,   0, 1, 1, 1, 1, 1},
  { 1,   1, 1, 0, 1, 1, 1},
  { 0,   0, 1, 0, 0, 1, 1},
  { 1,   0, 0, 0, 1, 1, 1},
  { 2,   0, 0, 1, 1, 0, 1},
  { 2,   1, 1, 1, 1, 0, 1},
  { 2,   1, 1, 0, 0, 0, 1},
  { 2,   0, 0, 0, 0, 0, 1}
};

int Isosurface::EdgeData[256][13]=  
{
  {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 8, 3, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 9, 0, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 8, 3, 1, 8, 1, 9,-1,-1,-1,-1,-1,-1,-1},
  {10, 1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 8, 3, 0, 1, 2,10,-1,-1,-1,-1,-1,-1,-1},
  { 9, 0, 2, 9, 2,10,-1,-1,-1,-1,-1,-1,-1},
  { 3, 2, 8, 2,10, 8, 8,10, 9,-1,-1,-1,-1},
  {11, 2, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  {11, 2, 0,11, 0, 8,-1,-1,-1,-1,-1,-1,-1},
  {11, 2, 3, 0, 1, 9,-1,-1,-1,-1,-1,-1,-1},
  { 2, 1,11, 1, 9,11,11, 9, 8,-1,-1,-1,-1},
  {10, 1, 3,10, 3,11,-1,-1,-1,-1,-1,-1,-1},
  { 1, 0,10, 0, 8,10,10, 8,11,-1,-1,-1,-1},
  { 0, 3, 9, 3,11, 9, 9,11,10,-1,-1,-1,-1},
  { 8,10, 9, 8,11,10,-1,-1,-1,-1,-1,-1,-1},
  { 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 3, 0, 4, 3, 4, 7,-1,-1,-1,-1,-1,-1,-1},
  { 1, 9, 0, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1},
  { 9, 4, 1, 4, 7, 1, 1, 7, 3,-1,-1,-1,-1},
  {10, 1, 2, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1},
  { 2,10, 1, 0, 4, 7, 0, 7, 3,-1,-1,-1,-1},
  { 4, 7, 8, 0, 2,10, 0,10, 9,-1,-1,-1,-1},
  { 2, 7, 3, 2, 9, 7, 7, 9, 4, 2,10, 9,-1},
  { 2, 3,11, 7, 8, 4,-1,-1,-1,-1,-1,-1,-1},
  { 7,11, 4,11, 2, 4, 4, 2, 0,-1,-1,-1,-1},
  { 3,11, 2, 4, 7, 8, 9, 0, 1,-1,-1,-1,-1},
  { 2, 7,11, 2, 1, 7, 1, 4, 7, 1, 9, 4,-1},
  { 8, 4, 7,11,10, 1,11, 1, 3,-1,-1,-1,-1},
  {11, 4, 7, 1, 4,11, 1,11,10, 1, 0, 4,-1},
  { 3, 8, 0, 7,11, 4,11, 9, 4,11,10, 9,-1},
  { 7,11, 4, 4,11, 9,11,10, 9,-1,-1,-1,-1},
  { 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 3, 0, 8, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1},
  { 5, 4, 0, 5, 0, 1,-1,-1,-1,-1,-1,-1,-1},
  { 4, 8, 5, 8, 3, 5, 5, 3, 1,-1,-1,-1,-1},
  { 2,10, 1, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1},
  { 0, 8, 3, 5, 4, 9,10, 1, 2,-1,-1,-1,-1},
  {10, 5, 2, 5, 4, 2, 2, 4, 0,-1,-1,-1,-1},
  { 3, 4, 8, 3, 2, 4, 2, 5, 4, 2,10, 5,-1},
  {11, 2, 3, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1},
  { 9, 5, 4, 8,11, 2, 8, 2, 0,-1,-1,-1,-1},
  { 3,11, 2, 1, 5, 4, 1, 4, 0,-1,-1,-1,-1},
  { 8, 5, 4, 2, 5, 8, 2, 8,11, 2, 1, 5,-1},
  { 5, 4, 9, 1, 3,11, 1,11,10,-1,-1,-1,-1},
  { 0, 9, 1, 4, 8, 5, 8,10, 5, 8,11,10,-1},
  { 3, 4, 0, 3,10, 4, 4,10, 5, 3,11,10,-1},
  { 4, 8, 5, 5, 8,10, 8,11,10,-1,-1,-1,-1},
  { 9, 5, 7, 9, 7, 8,-1,-1,-1,-1,-1,-1,-1},
  { 0, 9, 3, 9, 5, 3, 3, 5, 7,-1,-1,-1,-1},
  { 8, 0, 7, 0, 1, 7, 7, 1, 5,-1,-1,-1,-1},
  { 1, 7, 3, 1, 5, 7,-1,-1,-1,-1,-1,-1,-1},
  { 1, 2,10, 5, 7, 8, 5, 8, 9,-1,-1,-1,-1},
  { 9, 1, 0,10, 5, 2, 5, 3, 2, 5, 7, 3,-1},
  { 5, 2,10, 8, 2, 5, 8, 5, 7, 8, 0, 2,-1},
  {10, 5, 2, 2, 5, 3, 5, 7, 3,-1,-1,-1,-1},
  {11, 2, 3, 8, 9, 5, 8, 5, 7,-1,-1,-1,-1},
  { 9, 2, 0, 9, 7, 2, 2, 7,11, 9, 5, 7,-1},
  { 0, 3, 8, 2, 1,11, 1, 7,11, 1, 5, 7,-1},
  { 2, 1,11,11, 1, 7, 1, 5, 7,-1,-1,-1,-1},
  { 3, 9, 1, 3, 8, 9, 7,11,10, 7,10, 5,-1},
  { 9, 1, 0,10, 7,11,10, 5, 7,-1,-1,-1,-1},
  { 3, 8, 0, 7,10, 5, 7,11,10,-1,-1,-1,-1},
  {11, 5, 7,11,10, 5,-1,-1,-1,-1,-1,-1,-1},
  {10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 8, 3, 0,10, 6, 5,-1,-1,-1,-1,-1,-1,-1},
  { 0, 1, 9, 5,10, 6,-1,-1,-1,-1,-1,-1,-1},
  {10, 6, 5, 9, 8, 3, 9, 3, 1,-1,-1,-1,-1},
  { 1, 2, 6, 1, 6, 5,-1,-1,-1,-1,-1,-1,-1},
  { 0, 8, 3, 2, 6, 5, 2, 5, 1,-1,-1,-1,-1},
  { 5, 9, 6, 9, 0, 6, 6, 0, 2,-1,-1,-1,-1},
  { 9, 6, 5, 3, 6, 9, 3, 9, 8, 3, 2, 6,-1},
  { 3,11, 2,10, 6, 5,-1,-1,-1,-1,-1,-1,-1},
  { 6, 5,10, 2, 0, 8, 2, 8,11,-1,-1,-1,-1},
  { 1, 9, 0, 6, 5,10,11, 2, 3,-1,-1,-1,-1},
  { 1,10, 2, 5, 9, 6, 9,11, 6, 9, 8,11,-1},
  {11, 6, 3, 6, 5, 3, 3, 5, 1,-1,-1,-1,-1},
  { 0, 5, 1, 0,11, 5, 5,11, 6, 0, 8,11,-1},
  { 0, 5, 9, 0, 3, 5, 3, 6, 5, 3,11, 6,-1},
  { 5, 9, 6, 6, 9,11, 9, 8,11,-1,-1,-1,-1},
  {10, 6, 5, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1},
  { 5,10, 6, 7, 3, 0, 7, 0, 4,-1,-1,-1,-1},
  { 5,10, 6, 0, 1, 9, 8, 4, 7,-1,-1,-1,-1},
  { 4, 5, 9, 6, 7,10, 7, 1,10, 7, 3, 1,-1},
  { 7, 8, 4, 5, 1, 2, 5, 2, 6,-1,-1,-1,-1},
  { 4, 1, 0, 4, 5, 1, 6, 7, 3, 6, 3, 2,-1},
  { 9, 4, 5, 8, 0, 7, 0, 6, 7, 0, 2, 6,-1},
  { 4, 5, 9, 6, 3, 2, 6, 7, 3,-1,-1,-1,-1},
  { 7, 8, 4, 2, 3,11,10, 6, 5,-1,-1,-1,-1},
  {11, 6, 7,10, 2, 5, 2, 4, 5, 2, 0, 4,-1},
  {11, 6, 7, 8, 0, 3, 1,10, 2, 9, 4, 5,-1},
  { 6, 7,11, 1,10, 2, 9, 4, 5,-1,-1,-1,-1},
  { 6, 7,11, 4, 5, 8, 5, 3, 8, 5, 1, 3,-1},
  { 6, 7,11, 4, 1, 0, 4, 5, 1,-1,-1,-1,-1},
  { 4, 5, 9, 3, 8, 0,11, 6, 7,-1,-1,-1,-1},
  { 9, 4, 5, 7,11, 6,-1,-1,-1,-1,-1,-1,-1},
  {10, 6, 4,10, 4, 9,-1,-1,-1,-1,-1,-1,-1},
  { 8, 3, 0, 9,10, 6, 9, 6, 4,-1,-1,-1,-1},
  { 1,10, 0,10, 6, 0, 0, 6, 4,-1,-1,-1,-1},
  { 8, 6, 4, 8, 1, 6, 6, 1,10, 8, 3, 1,-1},
  { 9, 1, 4, 1, 2, 4, 4, 2, 6,-1,-1,-1,-1},
  { 1, 0, 9, 3, 2, 8, 2, 4, 8, 2, 6, 4,-1},
  { 2, 4, 0, 2, 6, 4,-1,-1,-1,-1,-1,-1,-1},
  { 3, 2, 8, 8, 2, 4, 2, 6, 4,-1,-1,-1,-1},
  { 2, 3,11, 6, 4, 9, 6, 9,10,-1,-1,-1,-1},
  { 0,10, 2, 0, 9,10, 4, 8,11, 4,11, 6,-1},
  {10, 2, 1,11, 6, 3, 6, 0, 3, 6, 4, 0,-1},
  {10, 2, 1,11, 4, 8,11, 6, 4,-1,-1,-1,-1},
  { 1, 4, 9,11, 4, 1,11, 1, 3,11, 6, 4,-1},
  { 0, 9, 1, 4,11, 6, 4, 8,11,-1,-1,-1,-1},
  {11, 6, 3, 3, 6, 0, 6, 4, 0,-1,-1,-1,-1},
  { 8, 6, 4, 8,11, 6,-1,-1,-1,-1,-1,-1,-1},
  { 6, 7,10, 7, 8,10,10, 8, 9,-1,-1,-1,-1},
  { 9, 3, 0, 6, 3, 9, 6, 9,10, 6, 7, 3,-1},
  { 6, 1,10, 6, 7, 1, 7, 0, 1, 7, 8, 0,-1},
  { 6, 7,10,10, 7, 1, 7, 3, 1,-1,-1,-1,-1},
  { 7, 2, 6, 7, 9, 2, 2, 9, 1, 7, 8, 9,-1},
  { 1, 0, 9, 3, 6, 7, 3, 2, 6,-1,-1,-1,-1},
  { 8, 0, 7, 7, 0, 6, 0, 2, 6,-1,-1,-1,-1},
  { 2, 7, 3, 2, 6, 7,-1,-1,-1,-1,-1,-1,-1},
  { 7,11, 6, 3, 8, 2, 8,10, 2, 8, 9,10,-1},
  {11, 6, 7,10, 0, 9,10, 2, 0,-1,-1,-1,-1},
  { 2, 1,10, 7,11, 6, 8, 0, 3,-1,-1,-1,-1},
  { 1,10, 2, 6, 7,11,-1,-1,-1,-1,-1,-1,-1},
  { 7,11, 6, 3, 9, 1, 3, 8, 9,-1,-1,-1,-1},
  { 9, 1, 0,11, 6, 7,-1,-1,-1,-1,-1,-1,-1},
  { 0, 3, 8,11, 6, 7,-1,-1,-1,-1,-1,-1,-1},
  {11, 6, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  {11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 0, 8, 3,11, 7, 6,-1,-1,-1,-1,-1,-1,-1},
  { 9, 0, 1,11, 7, 6,-1,-1,-1,-1,-1,-1,-1},
  { 7, 6,11, 3, 1, 9, 3, 9, 8,-1,-1,-1,-1},
  { 1, 2,10, 6,11, 7,-1,-1,-1,-1,-1,-1,-1},
  { 2,10, 1, 7, 6,11, 8, 3, 0,-1,-1,-1,-1},
  {11, 7, 6,10, 9, 0,10, 0, 2,-1,-1,-1,-1},
  { 7, 6,11, 3, 2, 8, 8, 2,10, 8,10, 9,-1},
  { 2, 3, 7, 2, 7, 6,-1,-1,-1,-1,-1,-1,-1},
  { 8, 7, 0, 7, 6, 0, 0, 6, 2,-1,-1,-1,-1},
  { 1, 9, 0, 3, 7, 6, 3, 6, 2,-1,-1,-1,-1},
  { 7, 6, 2, 7, 2, 9, 2, 1, 9, 7, 9, 8,-1},
  { 6,10, 7,10, 1, 7, 7, 1, 3,-1,-1,-1,-1},
  { 6,10, 1, 6, 1, 7, 7, 1, 0, 7, 0, 8,-1},
  { 9, 0, 3, 6, 9, 3, 6,10, 9, 6, 3, 7,-1},
  { 6,10, 7, 7,10, 8,10, 9, 8,-1,-1,-1,-1},
  { 8, 4, 6, 8, 6,11,-1,-1,-1,-1,-1,-1,-1},
  {11, 3, 6, 3, 0, 6, 6, 0, 4,-1,-1,-1,-1},
  { 0, 1, 9, 4, 6,11, 4,11, 8,-1,-1,-1,-1},
  { 1, 9, 4,11, 1, 4,11, 3, 1,11, 4, 6,-1},
  {10, 1, 2,11, 8, 4,11, 4, 6,-1,-1,-1,-1},
  {10, 1, 2,11, 3, 6, 6, 3, 0, 6, 0, 4,-1},
  { 0, 2,10, 0,10, 9, 4,11, 8, 4, 6,11,-1},
  { 2,11, 3, 6, 9, 4, 6,10, 9,-1,-1,-1,-1},
  { 3, 8, 2, 8, 4, 2, 2, 4, 6,-1,-1,-1,-1},
  { 2, 0, 4, 2, 4, 6,-1,-1,-1,-1,-1,-1,-1},
  { 1, 9, 0, 3, 8, 2, 2, 8, 4, 2, 4, 6,-1},
  { 9, 4, 1, 1, 4, 2, 4, 6, 2,-1,-1,-1,-1},
  { 8, 4, 6, 8, 6, 1, 6,10, 1, 8, 1, 3,-1},
  { 1, 0,10,10, 0, 6, 0, 4, 6,-1,-1,-1,-1},
  { 8, 0, 3, 9, 6,10, 9, 4, 6,-1,-1,-1,-1},
  {10, 4, 6,10, 9, 4,-1,-1,-1,-1,-1,-1,-1},
  { 9, 5, 4, 7, 6,11,-1,-1,-1,-1,-1,-1,-1},
  { 4, 9, 5, 3, 0, 8,11, 7, 6,-1,-1,-1,-1},
  { 6,11, 7, 4, 0, 1, 4, 1, 5,-1,-1,-1,-1},
  { 6,11, 7, 4, 8, 5, 5, 8, 3, 5, 3, 1,-1},
  { 6,11, 7, 1, 2,10, 9, 5, 4,-1,-1,-1,-1},
  {11, 7, 6, 8, 3, 0, 1, 2,10, 9, 5, 4,-1},
  {11, 7, 6,10, 5, 2, 2, 5, 4, 2, 4, 0,-1},
  { 7, 4, 8, 2,11, 3,10, 5, 6,-1,-1,-1,-1},
  { 4, 9, 5, 6, 2, 3, 6, 3, 7,-1,-1,-1,-1},
  { 9, 5, 4, 8, 7, 0, 0, 7, 6, 0, 6, 2,-1},
  { 4, 0, 1, 4, 1, 5, 6, 3, 7, 6, 2, 3,-1},
  { 7, 4, 8, 5, 2, 1, 5, 6, 2,-1,-1,-1,-1},
  { 4, 9, 5, 6,10, 7, 7,10, 1, 7, 1, 3,-1},
  { 5, 6,10, 0, 9, 1, 8, 7, 4,-1,-1,-1,-1},
  { 5, 6,10, 7, 0, 3, 7, 4, 0,-1,-1,-1,-1},
  {10, 5, 6, 4, 8, 7,-1,-1,-1,-1,-1,-1,-1},
  { 5, 6, 9, 6,11, 9, 9,11, 8,-1,-1,-1,-1},
  { 0, 9, 5, 0, 5, 3, 3, 5, 6, 3, 6,11,-1},
  { 0, 1, 5, 0, 5,11, 5, 6,11, 0,11, 8,-1},
  {11, 3, 6, 6, 3, 5, 3, 1, 5,-1,-1,-1,-1},
  { 1, 2,10, 5, 6, 9, 9, 6,11, 9,11, 8,-1},
  { 1, 0, 9, 6,10, 5,11, 3, 2,-1,-1,-1,-1},
  { 6,10, 5, 2, 8, 0, 2,11, 8,-1,-1,-1,-1},
  { 3, 2,11,10, 5, 6,-1,-1,-1,-1,-1,-1,-1},
  { 9, 5, 6, 3, 9, 6, 3, 8, 9, 3, 6, 2,-1},
  { 5, 6, 9, 9, 6, 0, 6, 2, 0,-1,-1,-1,-1},
  { 0, 3, 8, 2, 5, 6, 2, 1, 5,-1,-1,-1,-1},
  { 1, 6, 2, 1, 5, 6,-1,-1,-1,-1,-1,-1,-1},
  {10, 5, 6, 9, 3, 8, 9, 1, 3,-1,-1,-1,-1},
  { 0, 9, 1, 5, 6,10,-1,-1,-1,-1,-1,-1,-1},
  { 8, 0, 3,10, 5, 6,-1,-1,-1,-1,-1,-1,-1},
  {10, 5, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  {11, 7, 5,11, 5,10,-1,-1,-1,-1,-1,-1,-1},
  { 3, 0, 8, 7, 5,10, 7,10,11,-1,-1,-1,-1},
  { 9, 0, 1,10,11, 7,10, 7, 5,-1,-1,-1,-1},
  { 3, 1, 9, 3, 9, 8, 7,10,11, 7, 5,10,-1},
  { 2,11, 1,11, 7, 1, 1, 7, 5,-1,-1,-1,-1},
  { 0, 8, 3, 2,11, 1, 1,11, 7, 1, 7, 5,-1},
  { 9, 0, 2, 9, 2, 7, 2,11, 7, 9, 7, 5,-1},
  {11, 3, 2, 8, 5, 9, 8, 7, 5,-1,-1,-1,-1},
  {10, 2, 5, 2, 3, 5, 5, 3, 7,-1,-1,-1,-1},
  { 5,10, 2, 8, 5, 2, 8, 7, 5, 8, 2, 0,-1},
  { 9, 0, 1,10, 2, 5, 5, 2, 3, 5, 3, 7,-1},
  { 1,10, 2, 5, 8, 7, 5, 9, 8,-1,-1,-1,-1},
  { 1, 3, 7, 1, 7, 5,-1,-1,-1,-1,-1,-1,-1},
  { 8, 7, 0, 0, 7, 1, 7, 5, 1,-1,-1,-1,-1},
  { 0, 3, 9, 9, 3, 5, 3, 7, 5,-1,-1,-1,-1},
  { 9, 7, 5, 9, 8, 7,-1,-1,-1,-1,-1,-1,-1},
  { 4, 5, 8, 5,10, 8, 8,10,11,-1,-1,-1,-1},
  { 3, 0, 4, 3, 4,10, 4, 5,10, 3,10,11,-1},
  { 0, 1, 9, 4, 5, 8, 8, 5,10, 8,10,11,-1},
  { 5, 9, 4, 1,11, 3, 1,10,11,-1,-1,-1,-1},
  { 8, 4, 5, 2, 8, 5, 2,11, 8, 2, 5, 1,-1},
  { 3, 2,11, 1, 4, 5, 1, 0, 4,-1,-1,-1,-1},
  { 9, 4, 5, 8, 2,11, 8, 0, 2,-1,-1,-1,-1},
  {11, 3, 2, 9, 4, 5,-1,-1,-1,-1,-1,-1,-1},
  { 3, 8, 4, 3, 4, 2, 2, 4, 5, 2, 5,10,-1},
  {10, 2, 5, 5, 2, 4, 2, 0, 4,-1,-1,-1,-1},
  { 0, 3, 8, 5, 9, 4,10, 2, 1,-1,-1,-1,-1},
  { 2, 1,10, 9, 4, 5,-1,-1,-1,-1,-1,-1,-1},
  { 4, 5, 8, 8, 5, 3, 5, 1, 3,-1,-1,-1,-1},
  { 5, 0, 4, 5, 1, 0,-1,-1,-1,-1,-1,-1,-1},
  { 3, 8, 0, 4, 5, 9,-1,-1,-1,-1,-1,-1,-1},
  { 9, 4, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 7, 4,11, 4, 9,11,11, 9,10,-1,-1,-1,-1},
  { 3, 0, 8, 7, 4,11,11, 4, 9,11, 9,10,-1},
  {11, 7, 4, 1,11, 4, 1,10,11, 1, 4, 0,-1},
  { 8, 7, 4,11, 1,10,11, 3, 1,-1,-1,-1,-1},
  { 2,11, 7, 2, 7, 1, 1, 7, 4, 1, 4, 9,-1},
  { 3, 2,11, 4, 8, 7, 9, 1, 0,-1,-1,-1,-1},
  { 7, 4,11,11, 4, 2, 4, 0, 2,-1,-1,-1,-1},
  { 2,11, 3, 7, 4, 8,-1,-1,-1,-1,-1,-1,-1},
  { 2, 3, 7, 2, 7, 9, 7, 4, 9, 2, 9,10,-1},
  { 4, 8, 7, 0,10, 2, 0, 9,10,-1,-1,-1,-1},
  { 2, 1,10, 0, 7, 4, 0, 3, 7,-1,-1,-1,-1},
  {10, 2, 1, 8, 7, 4,-1,-1,-1,-1,-1,-1,-1},
  { 9, 1, 4, 4, 1, 7, 1, 3, 7,-1,-1,-1,-1},
  { 1, 0, 9, 8, 7, 4,-1,-1,-1,-1,-1,-1,-1},
  { 3, 4, 0, 3, 7, 4,-1,-1,-1,-1,-1,-1,-1},
  { 8, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 8, 9,10, 8,10,11,-1,-1,-1,-1,-1,-1,-1},
  { 0, 9, 3, 3, 9,11, 9,10,11,-1,-1,-1,-1},
  { 1,10, 0, 0,10, 8,10,11, 8,-1,-1,-1,-1},
  {10, 3, 1,10,11, 3,-1,-1,-1,-1,-1,-1,-1},
  { 2,11, 1, 1,11, 9,11, 8, 9,-1,-1,-1,-1},
  {11, 3, 2, 0, 9, 1,-1,-1,-1,-1,-1,-1,-1},
  {11, 0, 2,11, 8, 0,-1,-1,-1,-1,-1,-1,-1},
  {11, 3, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 3, 8, 2, 2, 8,10, 8, 9,10,-1,-1,-1,-1},
  { 9, 2, 0, 9,10, 2,-1,-1,-1,-1,-1,-1,-1},
  { 8, 0, 3, 1,10, 2,-1,-1,-1,-1,-1,-1,-1},
  {10, 2, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 8, 1, 3, 8, 9, 1,-1,-1,-1,-1,-1,-1,-1},
  { 9, 1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 8, 0, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}
};
