#include "hdrwriter.h"
#include <vector>
#include <string>
#include <assert.h>
#include <cmath>

typedef  unsigned int uint32_t;

// http://paulbourke.net/dataformats/tiff/ and http://partners.adobe.com/public/developer/en/tiff/TIFF6.pdf
// Not all programs support floating-point TIFFs, this was tested reading it back using Picturenaut and HDRShop
std::string EncodeFloatTIFF(unsigned int w, unsigned int h, float* RGBdata, unsigned int floatsPerPixel)
{
 assert(floatsPerPixel>=3); // we write only three floats (RGB) but support larger strides
 
 std::string outData;
 
 unsigned int image_size_bytes = w*h*3 * sizeof(float);
 outData.reserve(image_size_bytes + 500); // 500 is some slack for headers etc, I should compute it exactly... :)
 
 // Header
 outData.push_back(0x4d); outData.push_back(0x4d); // First two chars specify MM for big endian TODO - convert to little to make it easier on x86
 outData.push_back(0); outData.push_back(42); // Tiff version ID 
 
 unsigned int IFD_offset = 8 + image_size_bytes; // IFD table usually follows image
 outData.push_back((IFD_offset & 0xff000000) >> 24);
 outData.push_back((IFD_offset & 0xff0000) >> 16);
 outData.push_back((IFD_offset & 0xff00) >> 8);
 outData.push_back(IFD_offset & 0xff);
 
 // Image data
 for (unsigned int y=0; y<h; y++) 
 {
  for (unsigned int x=0; x<w; x++) 
  {
   unsigned int f = 0;
   for(; f< 3; f++,RGBdata++)
   {
    uint32_t floatAsInt = *reinterpret_cast<uint32_t*>(RGBdata);
    outData.push_back((floatAsInt & 0xff000000) >> 24);
    outData.push_back((floatAsInt & 0xff0000) >> 16);
    outData.push_back((floatAsInt & 0xff00) >> 8);
    outData.push_back(floatAsInt & 0xff);
   }
   for(; f<floatsPerPixel; f++)
    RGBdata++;
  }
 }
 
 // IFD Tags
 unsigned int NUM_IFD = 12;
 
 assert(outData.size() == IFD_offset);
 outData.push_back(0);
 outData.push_back(NUM_IFD); // Number of tags
 
 outData.push_back(1); outData.push_back(0); // -- width tag
 outData.push_back(0); outData.push_back(3); // short format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(1); // single value
 outData.push_back((w & 0xff00) >> 8); outData.push_back(w & 0xff); // value
 outData.push_back(0); outData.push_back(0); // padding (as we specified short value)
 
 outData.push_back(1); outData.push_back(1); // -- height tag
 outData.push_back(0); outData.push_back(3); // short format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(1); // single value
 outData.push_back((h & 0xff00) >> 8); outData.push_back(h & 0xff); // value
 outData.push_back(0); outData.push_back(0); // padding (as we specified short value)
 
 outData.push_back(1); outData.push_back(3); // -- compression
 outData.push_back(0); outData.push_back(3); // short format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(1); // single value
 outData.push_back(0); outData.push_back(1); // none
 outData.push_back(0); outData.push_back(0); // padding (as we specified short value)
 
 outData.push_back(1); outData.push_back(6); // -- photometric interpretation
 outData.push_back(0); outData.push_back(3); // short format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(1); // single value
 outData.push_back(0); outData.push_back(2); // RGB
 outData.push_back(0); outData.push_back(0); // padding (as we specified short value)
 
 outData.push_back(1); outData.push_back(0x12); // -- orientation
 outData.push_back(0); outData.push_back(3); // short format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(1); // single value
 outData.push_back(0); outData.push_back(1); // 
 outData.push_back(0); outData.push_back(0); // padding (as we specified short value)
 
 outData.push_back(1); outData.push_back(0x15); // -- samples per pixel
 outData.push_back(0); outData.push_back(3); // short format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(1); // single value
 outData.push_back(0); outData.push_back(3); // three samples (RGB)
 outData.push_back(0); outData.push_back(0); // padding (as we specified short value)
 
 outData.push_back(1); outData.push_back(0x16); // -- rows per strip
 outData.push_back(0); outData.push_back(3); // short format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(1); // single value
 outData.push_back((h & 0xff00) >> 8); outData.push_back(h & 0xff); // value
 outData.push_back(0); outData.push_back(0); // padding (as we specified short value)
 
 outData.push_back(1); outData.push_back(0x17); // -- strip byte count (total size)
 outData.push_back(0); outData.push_back(4); // long format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(1); // single value
 outData.push_back((image_size_bytes & 0xff000000) >> 24);
 outData.push_back((image_size_bytes & 0xff0000) >> 16);
 outData.push_back((image_size_bytes & 0xff00) >> 8);
 outData.push_back(image_size_bytes & 0xff);
 
 outData.push_back(1); outData.push_back(0x1c); // -- planar configuration
 outData.push_back(0); outData.push_back(3); // short format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(1); // single value
 outData.push_back(0); outData.push_back(1); // single image plane
 outData.push_back(0); outData.push_back(0); // padding (as we specified short value)
 
 outData.push_back(1); outData.push_back(0x11); // -- strip offset
 outData.push_back(0); outData.push_back(4); // long format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(1); // single value
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(8); // image starts right after the 8-byte header
 
 outData.push_back(1); outData.push_back(2); // -- bits per sample
 outData.push_back(0); outData.push_back(3); // short format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(3); // three values
 unsigned int BPS_offset = 8 + image_size_bytes + 2 + (NUM_IFD * 12) + 4; // offset to data (as data is > 4 bytes)
 outData.push_back((BPS_offset & 0xff000000) >> 24);
 outData.push_back((BPS_offset & 0xff0000) >> 16);
 outData.push_back((BPS_offset & 0xff00) >> 8);
 outData.push_back(BPS_offset & 0xff);
 
 outData.push_back(1); outData.push_back(0x53); // -- sample format
 outData.push_back(0); outData.push_back(3); // short format
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(3); // three values
 unsigned int SF_offset = BPS_offset + 3*2; // offset to data (as data is > 4 bytes)
 outData.push_back((SF_offset & 0xff000000) >> 24);
 outData.push_back((SF_offset & 0xff0000) >> 16);
 outData.push_back((SF_offset & 0xff00) >> 8);
 outData.push_back(SF_offset & 0xff);
 
 outData.push_back(0); outData.push_back(0); outData.push_back(0); outData.push_back(0); // IFD END
 
 // bits per sample data
 assert(outData.size() == BPS_offset);
 outData.push_back(0); outData.push_back(8*sizeof(float)); outData.push_back(0); outData.push_back(8*sizeof(float)); outData.push_back(0); outData.push_back(8*sizeof(float));
 
 // sample format data (1 = uint, 2 = sint, 3 = float)
 assert(outData.size() == SF_offset);
 outData.push_back(0); outData.push_back(3); outData.push_back(0); outData.push_back(3); outData.push_back(0); outData.push_back(3);
 
 return outData;
}

std::string EncodeRadianceHDR(unsigned int w, unsigned int h, const float* RGBdata, unsigned int floatsPerPixel)
{
assert(floatsPerPixel >= 3); // we write only three floats (RGB) but support larger strides

// Key-Value pairs after RADIANCE are optional
 //const char header[] = "#?RADIANCE\nEXPOSURE=1\nGAMMA=2.2\nFORMAT=32-bit_rle_rgbe\n\n";
 const char header[] = "#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n";
 
 std::string outData;
 outData.reserve(w*h*4 + sizeof(header) + 200); // 200 is some slack...
 
 std::vector<unsigned char> scanline[4];
 scanline[0].resize(w); scanline[1].resize(w); scanline[2].resize(w); scanline[3].resize(w);
 
 outData.append(header, sizeof(header)-1);
 outData.append("-Y ", 3); outData += std::to_string(h); outData.append(" +X ", 4); outData += std::to_string(w);
 outData.push_back('\n');
 
 for(unsigned int y=0 ; y<h; y++)
 {
  // RLE header
  // TODO looking at stb_image there seems to be also a non RLE line mode, which we should use as we don't really encode RLE here,
  // but I'm not sure that the way stb_image decodes the line header is standard-compliant...
  outData.push_back(2); outData.push_back(2);
  outData.push_back( (unsigned char)((w>>8) & 0xff) );
  outData.push_back( (unsigned char)(w & 0xff) );
 
  for(unsigned int x=0 ; x<w; x++)
  {
   unsigned char encodedPixel[4];
   float r = RGBdata[0], g = RGBdata[1], b= RGBdata[2];
 
   //r /= 179.0; g /= 179.0; b /= 179.0;  
   double maxV = r;
   if(maxV < g) maxV = g;
   if(maxV < b) maxV = b;
 
   if(maxV < std::numeric_limits<double>::epsilon())
   {
    encodedPixel[0] = encodedPixel[1] = encodedPixel[2] = encodedPixel[3] = 0;
   }
   else
   {
    int e;
    maxV = frexp(maxV, &e) * 256.0/maxV;
    encodedPixel[0] = (unsigned char)(maxV * r);
    encodedPixel[1] = (unsigned char)(maxV * g);
    encodedPixel[2] = (unsigned char)(maxV * b);
    encodedPixel[3] = (unsigned char)(e + 128);
   }
 
   scanline[0][x] = encodedPixel[0];
   scanline[1][x] = encodedPixel[1];
   scanline[2][x] = encodedPixel[2];
   scanline[3][x] = encodedPixel[3];
   RGBdata += floatsPerPixel;
  }
 
  // For simplicity, write all as it was not RLE...
  for(unsigned int line=0; line < 4; line++)
  {
   auto scanIter = scanline[line].begin();
   auto scanEnd = scanline[line].end();
  
   while( scanIter < scanEnd )
   {
    size_t remaining = scanEnd-scanIter;
    // the last bit in a char, if set, would indicate a RLE run, we want to avoid that
    unsigned char toWrite = remaining>127 ? 127 : (unsigned char)remaining; 
 
    outData.push_back(toWrite); // length of the "non run" data
    outData.append((char*)& scanIter[0], (size_t)toWrite);
 
    scanIter += (size_t)toWrite;
   }
  }
 }
 
 return outData;
}
