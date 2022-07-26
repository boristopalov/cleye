#include <stdio.h>
#include <stdlib.h>
#include <jpeglib.h>
#include <png.h>
#include <sys/ioctl.h>
#include <signal.h>
#include <jerror.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>


struct winsize ws;
struct img1d {
	unsigned int numComponents;
	unsigned int width, height;
	unsigned char* pixels;
};

struct img2d {
	unsigned int rows;
	unsigned int cols;
	unsigned int numComponents;
	unsigned char** pixels;
};

unsigned char xterm256_rgb_table[240][3] =  { 
  {0, 0, 0},
  {0, 0, 95},
  {0, 0, 135},
  {0, 0, 175},
  {0, 0, 215},
  {0, 0, 255},
  {0, 95, 0},
  {0, 95, 95},
  {0, 95, 135},
  {0, 95, 175},
  {0, 95, 215},
  {0, 95, 255},
  {0, 135, 0},
  {0, 135, 95},
  {0, 135, 135},
  {0, 135, 175},
  {0, 135, 215},
  {0, 135, 255},
  {0, 175, 0},
  {0, 175, 95},
  {0, 175, 135},
  {0, 175, 175},
  {0, 175, 215},
  {0, 175, 255},
  {0, 215, 0},
  {0, 215, 95},
  {0, 215, 135},
  {0, 215, 175},
  {0, 215, 215},
  {0, 215, 255},
  {0, 255, 0},
  {0, 255, 95},
  {0, 255, 135},
  {0, 255, 175},
  {0, 255, 215},
  {0, 255, 255},
  {95, 0, 0},
  {95, 0, 95},
  {95, 0, 135},
  {95, 0, 175},
  {95, 0, 215},
  {95, 0, 255},
  {95, 95, 0},
  {95, 95, 95},
  {95, 95, 135},
  {95, 95, 175},
  {95, 95, 215},
  {95, 95, 255},
  {95, 135, 0},
  {95, 135, 95},
  {95, 135, 135},
  {95, 135, 175},
  {95, 135, 215},
  {95, 135, 255},
  {95, 175, 0},
  {95, 175, 95},
  {95, 175, 135},
  {95, 175, 175},
  {95, 175, 215},
  {95, 175, 255},
  {95, 215, 0},
  {95, 215, 95},
  {95, 215, 135},
  {95, 215, 175},
  {95, 215, 215},
  {95, 215, 255},
  {95, 255, 0},
  {95, 255, 95},
  {95, 255, 135},
  {95, 255, 175},
  {95, 255, 215},
  {95, 255, 255},
  {135, 0, 0},
  {135, 0, 95},
  {135, 0, 135},
  {135, 0, 175},
  {135, 0, 215},
  {135, 0, 255},
  {135, 95, 0},
  {135, 95, 95},
  {135, 95, 135},
  {135, 95, 175},
  {135, 95, 215},
  {135, 95, 255},
  {135, 135, 0},
  {135, 135, 95},
  {135, 135, 135},
  {135, 135, 175},
  {135, 135, 215},
  {135, 135, 255},
  {135, 175, 0},
  {135, 175, 95},
  {135, 175, 135},
  {135, 175, 175},
  {135, 175, 215},
  {135, 175, 255},
  {135, 215, 0},
  {135, 215, 95},
  {135, 215, 135},
  {135, 215, 175},
  {135, 215, 215},
  {135, 215, 255},
  {135, 255, 0},
  {135, 255, 95},
  {135, 255, 135},
  {135, 255, 175},
  {135, 255, 215},
  {135, 255, 255},
  {175, 0, 0},
  {175, 0, 95},
  {175, 0, 135},
  {175, 0, 175},
  {175, 0, 215},
  {175, 0, 255},
  {175, 95, 0},
  {175, 95, 95},
  {175, 95, 135},
  {175, 95, 175},
  {175, 95, 215},
  {175, 95, 255},
  {175, 135, 0},
  {175, 135, 95},
  {175, 135, 135},
  {175, 135, 175},
  {175, 135, 215},
  {175, 135, 255},
  {175, 175, 0},
  {175, 175, 95},
  {175, 175, 135},
  {175, 175, 175},
  {175, 175, 215},
  {175, 175, 255},
  {175, 215, 0},
  {175, 215, 95},
  {175, 215, 135},
  {175, 215, 175},
  {175, 215, 215},
  {175, 215, 255},
  {175, 255, 0},
  {175, 255, 95},
  {175, 255, 135},
  {175, 255, 175},
  {175, 255, 215},
  {175, 255, 255},
  {215, 0, 0},
  {215, 0, 95},
  {215, 0, 135},
  {215, 0, 175},
  {215, 0, 215},
  {215, 0, 255},
  {215, 95, 0},
  {215, 95, 95},
  {215, 95, 135},
  {215, 95, 175},
  {215, 95, 215},
  {215, 95, 255},
  {215, 135, 0},
  {215, 135, 95},
  {215, 135, 135},
  {215, 135, 175},
  {215, 135, 215},
  {215, 135, 255},
  {215, 175, 0},
  {215, 175, 95},
  {215, 175, 135},
  {215, 175, 175},
  {215, 175, 215},
  {215, 175, 255},
  {215, 215, 0},
  {215, 215, 95},
  {215, 215, 135},
  {215, 215, 175},
  {215, 215, 215},
  {215, 215, 255},
  {215, 255, 0},
  {215, 255, 95},
  {215, 255, 135},
  {215, 255, 175},
  {215, 255, 215},
  {215, 255, 255},
  {255, 0, 0},
  {255, 0, 95},
  {255, 0, 135},
  {255, 0, 175},
  {255, 0, 215},
  {255, 0, 255},
  {255, 95, 0},
  {255, 95, 95},
  {255, 95, 135},
  {255, 95, 175},
  {255, 95, 215},
  {255, 95, 255},
  {255, 135, 0},
  {255, 135, 95},
  {255, 135, 135},
  {255, 135, 175},
  {255, 135, 215},
  {255, 135, 255},
  {255, 175, 0},
  {255, 175, 95},
  {255, 175, 135},
  {255, 175, 175},
  {255, 175, 215},
  {255, 175, 255},
  {255, 215, 0},
  {255, 215, 95},
  {255, 215, 135},
  {255, 215, 175},
  {255, 215, 215},
  {255, 215, 255},
  {255, 255, 0},
  {255, 255, 95},
  {255, 255, 135},
  {255, 255, 175},
  {255, 255, 215},
  {255, 255, 255},
  {8, 8, 8},
  {18, 18, 18},
  {28, 28, 28},
  {38, 38, 38},
  {48, 48, 48},
  {58, 58, 58},
  {68, 68, 68},
  {78, 78, 78},
  {88, 88, 88},
	{98, 98, 98}, 
  {108, 108, 108},
  {118, 118, 118},
  {128, 128, 128},
  {138, 138, 138},
  {148, 148, 148},
  {158, 158, 158},
  {168, 168, 168},
  {178, 178, 178},
  {188, 188, 188},
  {198, 198, 198},
  {208, 208, 208},
  {218, 218, 218},
  {228, 228, 228},
  {238, 238, 238}
};


int manhattanDistance(int x1, int y1, int z1, int x2, int y2, int z2) { 
  int xdist = abs(x2 - x1);
  int ydist = abs(y2 - y1);
  int zdist = abs(z2 - z1);
  return xdist + ydist + zdist;
}


struct img1d* loadImage(char* filename){
  struct jpeg_decompress_struct info;
  struct jpeg_error_mgr error;
  struct img1d* img = (struct img1d*) malloc(sizeof(struct img1d));
  
  FILE* f;
  f = fopen(filename , "r");
  if (!f) {
    printf("error opening %s\n", filename);
    exit(1);
  }

  info.err = jpeg_std_error(&error);
  jpeg_create_decompress(&info);
  jpeg_stdio_src(&info, f);

  (void) jpeg_read_header(&info, TRUE);
  (void) jpeg_start_decompress(&info);

  img->numComponents = info.num_components;
  img->width = info.output_width;
  img->height = info.output_height;
  unsigned int row_stride = info.output_width * info.num_components;
  img->pixels = (unsigned char*) malloc(row_stride * info.output_height);

  /* 
    * JSAMPARRAY: a pointer to rows (2d array)
    * JSAMPROW: a pointer to a row of pixel samples
    * JSAMPLE: unsigned char - 0...255
    
    * From https://github.com/LuaDist/libjpeg/blob/master/example.c:
    	* The standard input image format is a rectangular array of pixels, with
 	* each pixel having the same number of "component" values (color channels).
 	* Each pixel row is an array of JSAMPLEs (which typically are unsigned chars).
 	* If you are working with color data, then the color values for each pixel
 	* must be adjacent in the row; for example, R,G,B,R,G,B,R,G,B,... for 24-bit
 	* RGB color.
   */

  
  JSAMPARRAY buffer = (JSAMPARRAY) malloc(sizeof(unsigned char*) * 1);
  buffer[0] = (unsigned char*) malloc(sizeof(unsigned char) * row_stride);

  long num_pixels = 0;
  while (info.output_scanline < info.output_height) {
    (void) jpeg_read_scanlines(&info, buffer, 1);
    memcpy(img->pixels + num_pixels, buffer[0], row_stride);
    num_pixels += row_stride;

  }

  jpeg_finish_decompress(&info);
  jpeg_destroy_decompress(&info);
  fclose(f);
  return img;
}

void write_image(char* filename, struct img1d* img) {
	struct jpeg_compress_struct info;
	struct jpeg_error_mgr err;
	FILE* outfile;

	info.err = jpeg_std_error(&err);
	
	jpeg_create_compress(&info);
	if ((outfile = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    exit(1);
  }
	jpeg_stdio_dest(&info, outfile);

	info.image_width = img->width;
	info.image_height = img->height;
	info.input_components = img->numComponents;
	if (info.input_components == 3) { 
		info.in_color_space = JCS_RGB;
	}
	if (info.input_components == 1) { 
		info.in_color_space = JCS_GRAYSCALE;
	}

	jpeg_set_defaults(&info);

	jpeg_start_compress(&info, true);
	int row_stride = info.image_width * info.input_components;
  JSAMPARRAY buffer = (JSAMPARRAY) malloc(sizeof(unsigned char*) * 1);
  buffer[0] = (unsigned char*) malloc(sizeof(unsigned char) * row_stride);
	long num_pixels = 0;
	while (info.next_scanline < info.image_height) {
		memcpy(buffer[0], img->pixels + num_pixels, row_stride);
		(void) jpeg_write_scanlines(&info, buffer, 1);
		num_pixels += row_stride;
	}
	jpeg_finish_compress(&info);
	fclose(outfile);
	jpeg_destroy_compress(&info);
	printf("wrote image to %s\n", filename);
}


unsigned char xtermCode(int r, int g, int b) {
  int min_dist = 765; // max manhattan distance between (0,0,0) and (255, 255, 255)
  unsigned char closest_color_match_idx = -1;
  int xterm_r, xterm_g, xterm_b;

  for (int i = 0; i < 240; i++) {
    xterm_r = xterm256_rgb_table[i][0];
    xterm_g = xterm256_rgb_table[i][1];
    xterm_b = xterm256_rgb_table[i][2];
    float dist = manhattanDistance(r, g, b, xterm_r, xterm_g, xterm_b);
    if (dist < min_dist) { 
      min_dist = dist;
      closest_color_match_idx = i;
    }  
  }
  return closest_color_match_idx+16; //+16 to avoid intersecting colors with ANSI default terminal colors
}

void draw_image(struct img1d* img) {
  ioctl(0, TIOCGWINSZ, &ws);
	float scale;
	int width = ws.ws_col;
	int height = ws.ws_row;
	float terminal_aspect_ratio = width / height;
	float image_aspect_ratio = (float)img->width / (float)img->height;
	if (terminal_aspect_ratio > image_aspect_ratio) {
    scale = (float)height / (float)img->height;
    width = (int)((float)img->width * scale) * 2;
  } else {
    scale = (float)width / (float)img->width;
    height = (int)((float)img->height * scale) / 2;
  }

	if (img->height < height) { 
		height = img->height;
	}
	if (img->width < width) { 
		width = img->width;
	}

  for (int y = 0; y < height; y++) {
    int img_y_index = (int)(((float)y / (float)height) * (float)img->height); 
    int y_offset = img_y_index * img->width * img->numComponents;
    for (int x = 0; x < width; x++) { 
      int img_x_index = (int)(((float)x / (float)width) * (float)img->width);
      int x_offset = img_x_index * img->numComponents;
      int index = y_offset + x_offset;
      int r = img->pixels[index];
      int g = img->pixels[index+1];
      int b = img->pixels[index+2];
			unsigned char xterm_code;
			if (img->numComponents == 1) { 
				xterm_code = xtermCode(r, r, r);
			} 
			if (img->numComponents == 3) {
      	xterm_code = xtermCode(r, g, b);
			}
      printf("\033[48;5;%im ", xterm_code);
    }
   printf("\n");
  }
}

struct img1d* grayscaleFilter(struct img1d* img) {
  unsigned char r, g, b;

  for (int i=0; i < img->width*img->height*img->numComponents; i+=3) {
    r = img->pixels[i];
    g = img->pixels[i+1];
    b = img->pixels[i+2];
    unsigned char luma = (unsigned char)(
			0.299f * (float)r
			+ 0.587f * (float)g
			+ 0.114f * (float)b
		);
    img->pixels[i] = luma;
    img->pixels[i+1] = luma;
    img->pixels[i+2] = luma;
  }
    return img;
}


struct img1d* convert_to_grayscale(struct img1d* img) {
	if (img -> numComponents != 3) { 
		printf("You can only convert RGB images to grayscale\n");
		return NULL;
	}
	struct img1d* grayscale_img = (struct img1d*) malloc(sizeof(struct img1d));
	grayscale_img->numComponents = 1;
	grayscale_img->height = img->height;
	grayscale_img->width = img->width;
	grayscale_img->pixels = (unsigned char*) malloc(img->height * img->width * sizeof(unsigned char*));
  unsigned char r, g, b;
	int j = 0;
  for (int i=0; i < img->width*img->height*3; i+=3) {
    r = img->pixels[i];
    g = img->pixels[i+1];
    b = img->pixels[i+2];
    unsigned char luma = (unsigned char)(
			0.299f * (float)r
			+ 0.587f * (float)g
			+ 0.114f * (float)b
		);
    grayscale_img->pixels[j] = luma;
		j++;
  }
    return grayscale_img;
}

/* 
 * https://stackoverflow.com/questions/31631224/hacks-for-clamping-integer-to-0-255-and-doubles-to-0-0-1-0
 * Clamps the input to a 0 to 255 range.
 * @param v any int value
 * @return {@code v < 0 ? 0 : v > 255 ? 255 : v}
 */
unsigned char clampTo8Bit(int v) {
    // if out of range
    if ((v & ~0xFF) != 0) {
        // invert sign bit, shift to fill, then mask (generates 0 or 255)
        v = ((~v) >> 31) & 0xFF;
    }
    return (unsigned char)v;
}

struct img2d* prewitt_img(struct img2d* img) {
	struct img2d* output = (struct img2d*) malloc(sizeof(struct img2d));
	output->rows = img->rows - 2;
	output->cols = img->cols - 2;
	output->numComponents = img->numComponents;
	output->pixels = (unsigned char**) malloc(output->rows * sizeof(unsigned char*));
	for (int i=0; i < output->rows; i++) { 
		output->pixels[i] = (unsigned char* )malloc(output->cols);
	}
	
	// ignore the edges
	for (int i = 1; i < img->rows-1; i++) { 
		for (int j = 1; j < (img->cols * img->numComponents)-1; j++) { 
			unsigned char Xval = (unsigned char)clampTo8Bit((int)(img->pixels[i-1][j-1] + img->pixels[i-1][j] + img->pixels[i-1][j+1] - img->pixels[i+1][j-1] - img->pixels[i+1][j] - img->pixels[i+1][j+1]));
			unsigned char Yval = (unsigned char)clampTo8Bit((int)(img->pixels[i-1][j+1] + img->pixels[i][j+1] + img->pixels[i+1][j+1] - img->pixels[i-1][j-1] - img->pixels[i][j-1] - img->pixels[i+1][j-1]));
			unsigned char new_pixel_val = (unsigned char) clampTo8Bit((int)(sqrt(pow((double)Xval,2)+pow((double)Yval, 2))));
			output->pixels[i-1][j-1] = new_pixel_val;
		}
	}
	return output;
}

struct img2d* gaussian_blur(struct img2d* img, int radius, float sigma) { 
	if (radius < 1) { 
		printf("radius must be at least 1\n");
		return NULL;
	}
	struct img2d* output = (struct img2d*) malloc(sizeof(struct img2d));
	output->rows = img->rows - 2;
	output->cols = img->cols - 2;
	output->numComponents = img->numComponents;
	output->pixels = (unsigned char**) malloc(output->rows * sizeof(unsigned char*));
	for (int i=0; i < output->rows; i++) { 
		output->pixels[i] = (unsigned char* )malloc(output->cols * output->numComponents);
	}
	
	int kernel_width = 2*radius+1;
	double** kernel = (double**) calloc(kernel_width, sizeof(double*));
	for (int i=0; i < kernel_width; i++) { 
		kernel[i] = (double*) calloc(kernel_width, sizeof(double));
	}
	double kernel_sum = 0.0;
	for (int x=-radius; x<=radius; x++) { 
		for (int y=-radius; y<=radius; y++) { 
			double numerator = (double)(-1 * (x * x + y * y));
			double denom = (double)2 * sigma * sigma;
			double fraction = numerator / denom;
			double epow = exp(fraction);
			double kernel_value = epow / (2 * M_PI * sigma * sigma);
			kernel_sum += kernel_value;
			kernel[x+radius][y+radius] = kernel_value; 
		}
	}
	for (int x=0; x<kernel_width; x++) { 
		for (int y=0; y<kernel_width; y++) { 
			kernel[x][y] /= kernel_sum; // normalize
		}
	}
	for (int i = radius; i < img->rows-radius; i++) { 
		for (int j = radius; j < (img->cols-radius) * img->numComponents; j+=img->numComponents) {

			double r = 0, g = 0, b = 0;
			for (int kx=-radius; kx<=radius; kx++) { 
				for (int ky=-radius; ky<=radius; ky++) { 
					double kernel_val = kernel[kx+radius][ky+radius];
					int imgx = i+kx;
					int imgy = j+ky;
					r += img->pixels[imgx][imgy] * kernel_val;
					if (img->numComponents == 1) { 
						continue;
					}
					g += img->pixels[imgx][imgy+1] * kernel_val;
					b += img->pixels[imgx][imgy+2] * kernel_val;
					//printf("kx: %d\n", kx);
				}
			}
			output->pixels[i-1][j-1] = clampTo8Bit((int)r);
			if (img->numComponents == 1) { 
				continue;
			}
			output->pixels[i-1][j] = clampTo8Bit((int)g);
			output->pixels[i-1][j+1] = clampTo8Bit((int)b);
		}
	}
	return output;
}

struct img2d* img1dto2d(struct img1d* img) {
	int k=0;
	int rows = img->height; 
	int cols = img->width * img->numComponents;
	unsigned char** new_2d = (unsigned char** )malloc(rows * sizeof(unsigned char* ));
	for (int i=0; i < rows; i++) {
		new_2d[i] = (unsigned char* )malloc(cols * sizeof(unsigned char));
		for (int j=0; j < cols; j++) { 
			new_2d[i][j] = img->pixels[k];
			k++;
		}
	}
	struct img2d* ret = (struct img2d*) malloc(sizeof(struct img2d));
	ret->numComponents = img->numComponents;
	ret->pixels = new_2d;
	ret->rows = rows;
	ret->cols = img->width;
	return ret;
}


struct img1d* img2dto1d(struct img2d* img_2d) { 
  struct img1d* img_1d = (struct img1d*) malloc(sizeof(struct img1d));
	img_1d->numComponents = img_2d->numComponents;
	img_1d->height = img_2d->rows;
	img_1d->width = img_2d->cols;
	int k = 0;

	img_1d->pixels = (unsigned char*) malloc(img_1d->numComponents * img_1d->height * img_1d->width);
	for (int i=0; i < img_2d->rows; i++) { 
		for (int j=0; j < img_2d->cols * img_2d->numComponents; j++) {
			img_1d->pixels[k] = img_2d->pixels[i][j];
			k++;
		}
	}
	return img_1d;
} 


int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("invalid arguments");
    exit(1);
  }

  struct img1d* img = loadImage(argv[1]);
	struct img2d* img_2d = img1dto2d(img);
	struct img2d* blurred_img = gaussian_blur(img_2d, 4, 1.0f);
	struct img1d* blurred_1d = img2dto1d(blurred_img);
	//draw_image(blurred_1d);
	bool grayscale_flag = false;
	if (argc == 3 && strcmp(argv[2], "-g") == 0) {
		grayscale_flag = true;
		//struct img1d* grayscale_img = convert_to_grayscale(img);
		//if (grayscale_img) { 		
			struct img2d* prewitt_2d = prewitt_img(blurred_img);
			struct img1d* prewitt_1d = img2dto1d(prewitt_2d);
			draw_image(prewitt_1d);
			write_image("lenaedge2.jpeg", prewitt_1d);
		
		//}
	}

  //draw_image(img);
  // signal( SIGWINCH, redraw );
  return 0;

   
}


