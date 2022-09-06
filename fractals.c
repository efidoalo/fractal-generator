/*================================;
 *
 * File: fractals.c
 * Content: Source code for a program
 * that lets you construct different fractals
 * with variable characteristics.
 * 
 * compile and link: gcc `pkg-config --cflags gtk+-3.0` -o fractals fractals.c `pkg-config --libs gtk+-3.0` -lm ~/Documents/Containers/C/vector.o -lmpfr -lgmp
 * Designed for little endian machines
 * Date: 8/5/2022
 *
 * 
 *********************************/

#include <gtk/gtk.h>
#include <cairo.h>
#include <math.h>
#include "../Containers/C/vector.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <mpfr.h>
#include <glib-object.h>

struct point
{
	double x;
	double y;
};

struct fpoint
{
	mpfr_t x;
	mpfr_t y;
};

// function that returns 2^n, where n is 0 or greater
unsigned int pow2(unsigned int n)
{
	unsigned int res = 1;
	for (unsigned int i=0; i<n; ++i) {
		res *= 2;
	}
	return res;
}

// Function that returns the 0-31 index bit position 
// (the right most bit position being 0) of the first 1
// (from right to left) found ijn mask.
// ie is mask = 0010, this function will return 1
unsigned char get_bit_position_from_mask(unsigned int mask)
{
	for (int i=0; i<32; ++i) {
		unsigned int temp_mask = 1 << i;
		if (mask & temp_mask) {
			return i;
		}
	}
}

// draw a fractal canopy with initial section
// of size length pixels. the angle is given in radians
// and determines the angle of each iteration in relation
// to the previous branch. ratio is given as a double between 0 and 1
// and defines the length of branches in relation to the length of the 
// branches in the previous iteration. Depth determines the number of iterations
// to carry out, and is at least one
void 
draw_fractal_canopy(cairo_t *cr, 
		    GtkStyleContext *context,
		    struct point *init_point,
		    double init_length, 
	            double angle,
		    double ratio,
		    int depth)
{
	unsigned int mask = 0;
	// create the first set of branches corresponding to zero value mask
	cairo_move_to(cr, init_point->x, init_point->y);
	cairo_line_to(cr, init_point->x, (init_point->y) - init_length);
	double current_length = init_length;
	struct point p1, p2;
	p1.x = init_point->x;
	p1.y = init_point->y;
	p2.x = init_point->x;
	p2.y = (init_point->y) - current_length;
	//cairo_stroke(cr);
	double angle1 = 0.0;
	for (int i=0; i<(depth-1); ++i) {
		double new_length = current_length*ratio;
		double angle2 = angle1 + angle;
		angle2 = fmod(angle2, 2.0*M_PI);
		if ( (fabs((M_PI/2.0) - angle2) < 1e-12) || 
		     (((M_PI/2.0) - angle2) > 0.0) ) {
			if (fabs((M_PI/2.0) - angle2) < 1e-12) {
				cairo_rel_line_to(cr, new_length, 0);
			}
			else {
				double angle3 = (M_PI/2.0) - angle2;
				double dx = new_length*cos(angle3);
				double dy = new_length*sin(angle3);
				cairo_rel_line_to(cr, dx, -dy);	
			}
		}
		else if ( (fabs(M_PI - angle2)<1e-12) ||
		          ((M_PI - angle2) > 0.0) ) {
			if (fabs(M_PI - angle2)<1e-12) {
				cairo_rel_line_to(cr, 0, new_length);
			}
			else {
				double angle3 = angle2 - (M_PI/2.0);
				double dx = new_length*cos(angle3);
				double dy = new_length*sin(angle3);
				cairo_rel_line_to(cr, dx, dy);
			}
		}
		else if ( (fabs((M_PI*(3.0/2.0)) - angle2)<1e-12) ||
			  (((M_PI*(3.0/2.0)) - angle2) > 0.0) ) {
			if (fabs((M_PI*(3.0/2.0)) - angle2)<1e-12) {
				cairo_rel_line_to(cr, -new_length, 0);
			}	
			else {
				double angle3 = angle2 - M_PI;
				double dx = new_length*sin(angle3);
				double dy = new_length*cos(angle3);
				cairo_rel_line_to(cr, -dx, dy);
			}
		}
		else if ( (fabs((M_PI*2.0) - angle2)<1e-12) ||
                          (((M_PI*2.0) - angle2) > 0.0) ) {
                        if (fabs((M_PI*2.0) - angle2)<1e-12) {
                                cairo_rel_line_to(cr, 0, -new_length);
                        }
                        else {
                                double angle3 = angle2 - (M_PI*(3.0/2.0));
                                double dx = new_length*cos(angle3);
                                double dy = new_length*sin(angle3);
                                cairo_rel_line_to(cr, -dx, -dy);
                        }
                }
		angle1 += angle;
		current_length = new_length;		
		 
	}
	cairo_stroke(cr);
	unsigned int NoOfEndNodes = pow2(depth-1);
	struct point divergence_point;
	divergence_point.x = p2.x;
	divergence_point.y = p2.y;

	for (int i=1; i<NoOfEndNodes; ++i) {
		++mask;
		angle1 = 0;
		current_length = init_length;
		unsigned char bit_position = get_bit_position_from_mask(mask);
		p2.x = divergence_point.x;
		p2.y = divergence_point.y;
		for (int j=depth-2; j>=0; --j) {
			current_length *= ratio;
			p1.x = p2.x;
			p1.y = p2.y;
			unsigned int temp_mask = 1 << j;
			if (mask & temp_mask) {
				angle1 += (2.0*M_PI) - angle;
			}
			else {
				angle1 += angle;
			}
			double angle_mod360 = fmod(angle1, 2*M_PI);
			if (fabs(angle_mod360)<1e-12) {
				p2.y = p2.y - current_length;
			}
			else if ((M_PI/2.0)-angle_mod360>0) {
				p2.x += sin(angle_mod360)*current_length;
				p2.y -= cos(angle_mod360)*current_length;
			}		
			else if (fabs((M_PI/2.0)-angle_mod360)<1e-12) {
				p2.x += current_length;
			}
			else if ((M_PI -angle_mod360)>0) {
				double temp_angle = angle_mod360 - (M_PI/2.0);
				p2.x += cos(temp_angle)*current_length;
				p2.y += sin(temp_angle)*current_length;
			}
			else if ( fabs(M_PI - angle_mod360)<1e-12) {
				p2.y += current_length;
			}
			else if ( ((M_PI*1.5) - angle_mod360) > 0 ) {
				double temp_angle = angle_mod360 - M_PI;
				p2.x -= sin(temp_angle)*current_length;
				p2.y += cos(temp_angle)*current_length;
			}
			else if ( fabs((M_PI*1.5)-angle_mod360)<1e-12) {
				p2.x -= current_length;
			}
			else {
				double temp_angle = angle_mod360 - (M_PI*1.5);
				p2.x -= cos(temp_angle)*current_length;
				p2.y -= sin(temp_angle)*current_length;
			}
			if (j==bit_position) {
				cairo_move_to(cr, p1.x, p1.y);	
			}
			if (j<=bit_position) {
				cairo_line_to(cr, p2.x, p2.y);
			}	
		}		
	}
	cairo_stroke(cr);

}

void draw_koch_curve(cairo_t *cr, 
		     GtkStyleContext *context, 
		     struct point *lower_left, 
		     int init_length, 
		     int iteration_depth)
{
	struct vector *prev_points = vector_null_init((int)sizeof(struct point));
	// construct  initial three points
	
	vector_push_back(prev_points, lower_left);
	struct point lower_right;
	lower_right.x = lower_left->x + init_length;
	lower_right.y = lower_left->y;
	vector_push_back(prev_points, &lower_right);
	struct point upper_center;
	upper_center.x = (lower_left->x + lower_right.x)/2;
	upper_center.y = ((double)lower_left->y) - ((double)init_length)*sin(M_PI/3.0);
	vector_push_back(prev_points, &upper_center);
	vector_push_back(prev_points, lower_left);

	for (int i=1; i<iteration_depth; ++i) {
		struct vector *curr_points_iter = vector_null_init((int)sizeof(struct point));
		int sizeof_prev_points = vector_get_size(prev_points);
		vector_push_back(curr_points_iter, lower_left);
		for (int j=0; j<sizeof_prev_points-1; ++j) {
			struct point *p0 = (struct point *)vector_read(prev_points, j);
			struct point *p1 = (struct point *)vector_read(prev_points, j+1);
		       	struct point p_third;
		        p_third.x =  p0->x + ((1.0/3.0)*((double)(p1->x - p0->x)));
			p_third.y =  p0->y + ((1.0/3.0)*((double)(p1->y - p0->y)));
			struct point p_two_thirds;
		        p_two_thirds.x = p0->x + ((2.0/3.0)*((double)(p1->x - p0->x)));	
			p_two_thirds.y = p0->y + ((2.0/3.0)*((double)(p1->y - p0->y)));
			int p_middle_third_length = sqrt(pow((p_two_thirds.x - p_third.x),2.0) + pow((p_two_thirds.y - p_third.y), 2.0));
			struct point p_outer;
			double bisecting_line_length = ((double)p_middle_third_length)*sin(M_PI/3.0);
			double scalar = bisecting_line_length/(double)p_middle_third_length;
			
			if (p_two_thirds.x != p_third.x) {
				p_outer.x = (((double)p_third.x + (double)p_two_thirds.x)/2.0) + (p_third.y - p_two_thirds.y)*scalar;
				p_outer.y = (((double)p_third.y + (double)p_two_thirds.y)/2.0) + (p_two_thirds.x - p_third.x)*scalar;
			}
			else if (p_two_thirds.y < p_third.y) {
				p_outer.x = p_two_thirds.x + bisecting_line_length;
				p_outer.y = (int)(((double )(p_two_thirds.y + p_third.y))/2.0);
			}
			else if (p_two_thirds.y > p_third.y) {
				p_outer.x = p_two_thirds.x - bisecting_line_length;
				p_outer.y = (int)(((double )(p_two_thirds.y + p_third.y))/2.0);
			}
			vector_push_back(curr_points_iter, &p_third);
			vector_push_back(curr_points_iter, &p_outer);
			vector_push_back(curr_points_iter, &p_two_thirds);
			vector_push_back(curr_points_iter, p1);
		}
		vector_free(prev_points);
		prev_points = curr_points_iter;
	}
	
	cairo_move_to(cr, lower_left->x, lower_left->y);
	for (int i=1; i<vector_get_size(prev_points); ++i) {
		struct point *curr_point = (struct point *)vector_read(prev_points,i );
		cairo_line_to(cr, curr_point->x, curr_point->y);
	}
	cairo_stroke(cr);
	
}

struct koch_curve_state
{
	int iterations;
	unsigned char valid_options_provided;
};
struct koch_curve_state kcs;

struct canopy_state
{
	int trunk_pixel_length;
	double angle;
	double ratio;
	int depth;	
	unsigned char valid_options_provided;
};
struct canopy_state cs;

struct mandelbrot_set_state
{
	struct fpoint tl; // the top left point of a rectangular subset of the mandelbrot set
	struct fpoint br;  // bottom right point of a rectangular subset of the mandelbrot set
	mpfr_t threshold;  // cut off threshold, if absolute value of the point after iterations is greater than threshold, it is not in the mandelbrot set 
	int iterations;	   // number of iterations to perform on each point in the rectangle
	struct fpoint curr_point;
	struct fpoint temp_point;
	mpfr_t x_interval;
	mpfr_t y_interval;
	mpfr_t temp_val1; // temp vals are required during mandelbrot set inclusion clauclation
	mpfr_t temp_val2;
	mpfr_t temp_val3; // used in loop
	unsigned char point_selected; // determines when a point has been selected that defines the tl or br of a zoom in rectangle
	unsigned char valid_options_provided;
};
struct mandelbrot_set_state mss;

void draw_mandelbrot_set(cairo_t *cr,
                         GtkStyleContext *context, 
		         int width,
		         int height)
{
	cairo_surface_t *original_surface = cairo_get_target(cr);

	cairo_rectangle_int_t surface_rect;
	surface_rect.x = 0;
	surface_rect.y = 25;
	surface_rect.width = 600;
	surface_rect.height = 500;
	cairo_surface_t *surface = cairo_surface_map_to_image(original_surface, &surface_rect); 
	cairo_surface_type_t type = cairo_surface_get_type(surface); // image type
		
	int surface_width = cairo_image_surface_get_width(surface);
	int surface_height = cairo_image_surface_get_height(surface);
	cairo_surface_flush(surface);
	unsigned char *data = cairo_image_surface_get_data(surface);
	int stride = cairo_image_surface_get_stride(surface);
	mpfr_sub(mss.x_interval, mss.br.x, mss.tl.x, MPFR_RNDN);
	
	mpfr_div_ui(mss.x_interval, mss.x_interval, 599, MPFR_RNDN);	
	mpfr_sub(mss.y_interval, mss.tl.y, mss.br.y, MPFR_RNDN);

	mpfr_div_ui(mss.y_interval, mss.y_interval, 499, MPFR_RNDN);
	if (mpfr_cmp(mss.x_interval, mss.y_interval)>0) {
		mpfr_set(mss.y_interval, mss.x_interval, MPFR_RNDN); 
	}	
	else {
		mpfr_set(mss.x_interval, mss.y_interval, MPFR_RNDN);
	}
	mpfr_set(mss.curr_point.x, mss.tl.x, MPFR_RNDN);
	mpfr_set(mss.curr_point.y, mss.tl.y, MPFR_RNDN);
	mpfr_set(mss.temp_val3, mss.curr_point.y, MPFR_RNDN);
	for (int i=0; i<600; ++i) {
		if (i!=0) {
			mpfr_add(mss.curr_point.x, mss.curr_point.x, mss.x_interval, MPFR_RNDN);
		}
		mpfr_set(mss.curr_point.y, mss.temp_val3, MPFR_RNDN);
		for (int j=0; j<500; ++j) {
			if (j!=0) {
				mpfr_sub(mss.curr_point.y, mss.curr_point.y, mss.y_interval, MPFR_RNDN);
			}
			mpfr_set_d(mss.temp_point.x, 0.0, MPFR_RNDN);
			mpfr_set_d(mss.temp_point.y, 0.0, MPFR_RNDN);
			unsigned char point_in_mandelbrot_set = 1;
			for (int k=0; k<mss.iterations; ++k) {
				if (k==0) {
					mpfr_set(mss.temp_point.x, mss.curr_point.x, MPFR_RNDN);
					mpfr_set(mss.temp_point.y, mss.curr_point.y, MPFR_RNDN);
				}
				else {
				// use temp_val1 and temp_val2 to update temp_point. temp_point = temp_point^2 + curr_point	
					mpfr_mul(mss.temp_val1, mss.temp_point.x, mss.temp_point.x, MPFR_RNDN);
					mpfr_set(mss.temp_val2, mss.temp_point.x, MPFR_RNDN);
					mpfr_set(mss.temp_point.x, mss.temp_val1, MPFR_RNDN);
					mpfr_mul(mss.temp_val1, mss.temp_point.y, mss.temp_point.y, MPFR_RNDN);
					mpfr_sub(mss.temp_point.x, mss.temp_point.x, mss.temp_val1, MPFR_RNDN);
					mpfr_mul_d(mss.temp_val2, mss.temp_val2, 2.0, MPFR_RNDN);
					mpfr_mul(mss.temp_point.y, mss.temp_point.y, mss.temp_val2, MPFR_RNDN);
					mpfr_add(mss.temp_point.x, mss.temp_point.x, mss.curr_point.x, MPFR_RNDN);
					mpfr_add(mss.temp_point.y, mss.temp_point.y, mss.curr_point.y, MPFR_RNDN);
				}
				mpfr_mul(mss.temp_val1, mss.temp_point.x, mss.temp_point.x, MPFR_RNDN);
				mpfr_mul(mss.temp_val2, mss.temp_point.y, mss.temp_point.y, MPFR_RNDN);
				mpfr_add(mss.temp_val1, mss.temp_val1, mss.temp_val2, MPFR_RNDN);
				mpfr_sqrt(mss.temp_val1, mss.temp_val1, MPFR_RNDN);
				if (mpfr_cmp(mss.temp_val1, mss.threshold) > 0) {
					point_in_mandelbrot_set = 0;
					break;
				}
			}
			if (point_in_mandelbrot_set) {
				data[(j*stride)+(i*4)] = 0;
				data[(j*stride)+(i*4)+1] = 0;
				data[(j*stride)+(i*4)+2] = 0;
                                data[(j*stride)+(i*4)+3] = 0;
			}		
		}
	}
	
	cairo_surface_mark_dirty(surface);
	cairo_surface_unmap_image(original_surface, surface);
}

unsigned int fractal_id; // integer value specifying which fractal to draw.
			 // -1 = no fractal specified (initial value)
			 // 0 = tree canopy,
			 // 1 = koch curve,
			 // 2 = mandelbrot set (zoom)

gulong button_press_hid = 0;

// this function is used to redraw the mandelbrot set with new boundary coordinates (zoomed in)
gboolean
button_press_callback(GtkWidget *drawing_area,
                      GdkEvent *event,
                      void *user_data)
{
        GdkEventButton *button_press_event = ((GdkEventButton *)event);
        if (mss.point_selected == 0) {
		// first point being selected (top left)
		mpfr_set(mss.curr_point.x, mss.tl.x, MPFR_RNDN);
		mpfr_set(mss.curr_point.y, mss.tl.y, MPFR_RNDN);
		mpfr_mul_ui(mss.temp_val1, mss.x_interval, (unsigned int)(button_press_event->x), MPFR_RNDN);
		mpfr_add(mss.curr_point.x, mss.curr_point.x, mss.temp_val1, MPFR_RNDN);
		mpfr_mul_ui(mss.temp_val1, mss.y_interval, (unsigned int)(button_press_event->y), MPFR_RNDN);
		mpfr_sub(mss.curr_point.y, mss.curr_point.y, mss.temp_val1, MPFR_RNDN);
	        mss.point_selected = 1;	
	}		
	else {
		// 2nd point being selected (bottom right)
		mpfr_set(mss.temp_point.x, mss.tl.x, MPFR_RNDN);
                mpfr_set(mss.temp_point.y, mss.tl.y, MPFR_RNDN);
                mpfr_mul_ui(mss.temp_val1, mss.x_interval, (unsigned int)(button_press_event->x), MPFR_RNDN);
                mpfr_add(mss.temp_point.x, mss.temp_point.x, mss.temp_val1, MPFR_RNDN);
                mpfr_mul_ui(mss.temp_val1, mss.y_interval, (unsigned int)(button_press_event->y), MPFR_RNDN);
                mpfr_sub(mss.temp_point.y, mss.temp_point.y, mss.temp_val1, MPFR_RNDN);
		mpfr_set(mss.tl.x, mss.curr_point.x, MPFR_RNDN);
		mpfr_set(mss.tl.y, mss.curr_point.y, MPFR_RNDN);
		mpfr_set(mss.br.x, mss.temp_point.x, MPFR_RNDN);
		mpfr_set(mss.br.y, mss.temp_point.y, MPFR_RNDN);
		gtk_widget_hide(drawing_area);
		gtk_widget_show(drawing_area);
		mss.point_selected = 0;
	}
}

unsigned char button_press_on_drawing_area_is_active(GtkWidget *drawing_area)
{
	if (g_signal_handler_is_connected(drawing_area, button_press_hid)) {
		return 1;	
	}
	else {
		return 0;
	}
}

gboolean 
draw_callback(GtkWidget *widget,
	      cairo_t *cr,
	      gpointer user_data)
{
	
	GtkStyleContext *context = gtk_widget_get_style_context(widget);
	int width = gtk_widget_get_allocated_width(widget);
	int height = gtk_widget_get_allocated_height(widget);
	
	switch (fractal_id) {
		case -1:
		{	// no fractal specified - do nothing
			if (button_press_on_drawing_area_is_active(widget)) {
				g_signal_handler_disconnect(widget, button_press_hid);
			}
			break;
			
		}
		case 0:  // tree canopy
		{
			double init_length = 100;
			double angle = M_PI/11.0;
			double ratio = 0.75;
			int depth = 12;
			struct point init_point;
			init_point.x = 300;
			init_point.y = 500;
			if (button_press_on_drawing_area_is_active(widget)) {
                                g_signal_handler_disconnect(widget, button_press_hid);
                        }
			draw_fractal_canopy(cr,
					    context,
					    &init_point,
					    cs.trunk_pixel_length,
					    cs.angle,
					    cs.ratio,
					    cs.depth);
			break;
		}
		case 1: 
		{       // koch curve
			int curve_width = 300;
			int curve_height = 400;
			int init_length = 300;
			struct point lower_left;
			lower_left.x = 150;
			lower_left.y = 350;
			int iterations = kcs.iterations;
			if (button_press_on_drawing_area_is_active(widget)) {
                                g_signal_handler_disconnect(widget, button_press_hid);
                        }
			draw_koch_curve(cr, context, &lower_left, init_length, iterations);
			break;
		}
		case 2:
		{	// draw mandelbrot set
			draw_mandelbrot_set(cr, context, width, height);
			if (!g_signal_handler_is_connected(widget, button_press_hid)) {
				button_press_hid = g_signal_connect(G_OBJECT(widget), "button-press-event", 
						                    G_CALLBACK(button_press_callback), 
								    NULL);
			}
			break;

		}

	}
	return FALSE;
}

// this function is used to redraw the mandelbrot set with new boundary coordinates (zoomed in)

void koch_curve_menu_item_activate(GtkMenuItem *menuitem, void *user_data)
{
	GtkStack *options_stack = (GtkStack *)user_data;
	fractal_id = -1;
	GtkWidget *koch_curve_options_grid = gtk_stack_get_child_by_name(options_stack, "koch curve options");
	gtk_widget_show(koch_curve_options_grid);
        gtk_stack_set_visible_child(options_stack, koch_curve_options_grid);
}

void canopy_menu_item_activate(GtkMenuItem *menuitem, void *user_data)
{
	fractal_id = -1;
	GtkWidget *options_stack = (GtkWidget *)user_data;
        GtkWidget *canopy_options_grid = gtk_stack_get_child_by_name(options_stack, "canopy options");
        gtk_widget_show(canopy_options_grid);
        gtk_stack_set_visible_child(options_stack, canopy_options_grid);
}

void mandelbrot_menu_item_activate(GtkMenuItem *menuitem, void *user_data)
{
	fractal_id = -1;
	GtkWidget *options_stack = (GtkWidget *)user_data;
	GtkWidget *mandelbrot_options_grid = gtk_stack_get_child_by_name(options_stack, "mandelbrot options");
	gtk_widget_show(mandelbrot_options_grid);
        gtk_stack_set_visible_child(options_stack, mandelbrot_options_grid);
}

void check_mandelbrot_options(GtkButton *button, void *os)
{
	GtkStack *options_stack = (GtkStack *)os;

	GtkGrid *mandelbrot_options_grid = (GtkGrid *)gtk_stack_get_child_by_name(options_stack, "mandelbrot options");
	GtkEntry *iterations_specifier = gtk_grid_get_child_at(mandelbrot_options_grid, 0, 0);
	GtkEntry *significand_length = gtk_grid_get_child_at(mandelbrot_options_grid, 0, 1);
	GtkEntry *threshold_specifier = gtk_grid_get_child_at(mandelbrot_options_grid, 0, 2);
	char *iterations = gtk_entry_get_text(iterations_specifier);
	int iterations_option_val = atoi(iterations);
	char *significand_length_str = gtk_entry_get_text(significand_length);
	int significand_len_option_val = atoi(significand_length_str);
	char *threshold_str = gtk_entry_get_text(threshold_specifier);
	double threshold_val = atof(threshold_str);
	if ((iterations_option_val > 0) && (significand_len_option_val >0) &&
	    (threshold_val > 3.0)) {
		mss.valid_options_provided = 1;
		mss.iterations = iterations_option_val;
		mss.point_selected = 0;
		mpfr_inits2(significand_len_option_val,
			    mss.tl.x, mss.tl.y, mss.br.x, mss.br.y, mss.threshold, mss.curr_point.x, mss.curr_point.y, mss.temp_point.x, mss.temp_point.y,
		    	    mss.x_interval, mss.y_interval, mss.temp_val1, mss.temp_val2, mss.temp_val3, (mpfr_ptr)0);
    		mpfr_init_set_d( mss.tl.x, -2.0, MPFR_RNDN);
		mpfr_init_set_d( mss.tl.y, 1.0, MPFR_RNDN);
		mpfr_init_set_d( mss.br.x, 1.0, MPFR_RNDN);
		mpfr_init_set_d( mss.br.y, -1.0, MPFR_RNDN);
		mpfr_init_set_d(mss.threshold, threshold_val, MPFR_RNDN);
		mpfr_init_set_d(mss.curr_point.x, -2.0, MPFR_RNDN);
		mpfr_init_set_d(mss.curr_point.y, 1.0, MPFR_RNDN);
		mpfr_init_set_d(mss.temp_point.x, 0.0, MPFR_RNDN);
		mpfr_init_set_d(mss.temp_point.y, 0.0, MPFR_RNDN);
		mpfr_init_set_d(mss.x_interval, 3.0, MPFR_RNDN);
		mpfr_init_set_d(mss.y_interval, 2.0, MPFR_RNDN);
		mpfr_init_set_d(mss.temp_val1, 0.0, MPFR_RNDN);
		mpfr_init_set_d(mss.temp_val2, 0.0, MPFR_RNDN);
		mpfr_init_set_d(mss.temp_val3, 0.0, MPFR_RNDN);		
	}
	else {
		mss.valid_options_provided = 0;
	}
}

void start_mandelbrot_draw(GtkButton *button, void *da)
{
	if (mss.valid_options_provided) {
		GtkDrawingArea *drawing_area = (GtkDrawingArea *)da;
		fractal_id = 2;
		gtk_widget_hide(drawing_area);
		gtk_widget_show(drawing_area);
	}
}

void check_canopy_options(GtkButton *button, void *os)
{
	GtkStack *options_stack = (GtkStack *)os;
	GtkGrid *canopy_options_grid = gtk_stack_get_child_by_name(options_stack, "canopy options");
	GtkEntry *trunk_length = gtk_grid_get_child_at(canopy_options_grid, 0, 0);
	GtkWidget *anglegrid = gtk_grid_get_child_at(canopy_options_grid, 0, 1);
	GtkScaleButton *angle_rad = gtk_grid_get_child_at(anglegrid, 1, 0);
	GtkEntry *ratio = gtk_grid_get_child_at(canopy_options_grid, 0, 2);
	GtkEntry *depth = gtk_grid_get_child_at(canopy_options_grid, 0, 3);
	char *trunk_pixel_len_str = gtk_entry_get_text(trunk_length);
	int trunk_len = atoi(trunk_pixel_len_str);
	char *ratio_str = gtk_entry_get_text(ratio);
	double ratio_val = atof(ratio_str);
	char *depth_str = gtk_entry_get_text(depth);
	int depth_val = atoi(depth_str);
	if (((0.0 < ratio_val) && (ratio_val < 1.0)) &&
	    (depth_val > 0) &&
	    (trunk_len > 0)) {	   
		cs.trunk_pixel_length = trunk_len;
		cs.angle = gtk_scale_button_get_value(angle_rad);
		cs.ratio = ratio_val;
		cs.depth = depth_val;
		cs.valid_options_provided = 1;
	}
	else {
		cs.valid_options_provided = 0;
	}
}

void start_draw_canopy(GtkButton *button, void *da)
{
	if (cs.valid_options_provided == 1) {
		fractal_id = 0; 
		GtkDrawingArea *drawing_area = (GtkDrawingArea *)da;
		gtk_widget_queue_draw(drawing_area);
	}
}

void check_koch_curve_options(GtkButton *button, void *os)
{
	GtkStack *options_stack = (GtkStack *)os;
	GtkGrid *koch_curve_options_grid = gtk_stack_get_child_by_name(options_stack, "koch curve options");
	GtkEntry *iterations_option = gtk_grid_get_child_at(koch_curve_options_grid, 0, 0);
	char *iterations_str = gtk_entry_get_text(iterations_option);
	int iteration_number = atoi(iterations_str);
	if (iteration_number > 0) {
		kcs.iterations = iteration_number;
		kcs.valid_options_provided = 1;	
	}
	else {
		kcs.valid_options_provided = 0;
	}
}

void start_koch_curve_draw(GtkButton *button, void *da)
{
	GtkDrawingArea *drawing_area = (GtkDrawingArea *)da;
	if (kcs.valid_options_provided == 1) {
		GtkDrawingArea *drawing_area = (GtkDrawingArea *)da;
		fractal_id = 1;
		gtk_widget_queue_draw(drawing_area);
	}
}

static void
activate (GtkApplication* app,
          gpointer        user_data)
{
  	GtkWidget *window;
	GtkWidget *main_grid;
	GtkWidget *menu_bar;
	GtkWidget *drawing_area;
	GtkWidget *fractal_selector_menu;
	GtkWidget *fractal_selector_menuitem;;
	GtkWidget *canopy_menu_item;
	GtkWidget *koch_curve_menu_item;
	GtkWidget *mandelbrot_menu_item;

	GtkWidget *options_stack;
	GtkWidget *options_texttagtable;
	GtkWidget *options_static_texttag; // not editable
	
	GtkWidget *mandelbrot_options_grid;
	GtkWidget *iterations_specifier;
	GtkWidget *significand_length_specifier;
	GtkWidget *threshold_specifier;
	GtkWidget *mandelbrot_start_button;
	GtkWidget *mandelbrot_text_buffer;
	GtkWidget *mandelbrot_text_view;

	GtkWidget *tree_canopy_options_grid;
	GtkWidget *initial_branch_length_specifier;
	GtkWidget *angle_specifier;
	GtkWidget *ratio;
	GtkWidget *depth;
	GtkWidget *canopy_start_button;
	GtkWidget *canopy_text_buffer;
	GtkWidget *canopy_text_view;	

	GtkWidget *koch_curve_options_grid;
	GtkWidget *koch_iterations_specifier;
	GtkWidget *koch_start_button;
	GtkWidget *koch_curve_text_buffer;
	GtkWidget *koch_curve_text_view;

	window = gtk_application_window_new (app);
  	gtk_window_set_title (GTK_WINDOW (window), "Fractal Generator");
  	gtk_window_set_default_size (GTK_WINDOW (window), 800, 600);
	
  	main_grid = gtk_grid_new();

	menu_bar = gtk_menu_bar_new();
	gtk_menu_bar_set_child_pack_direction(menu_bar, GTK_PACK_DIRECTION_TTB);
	options_stack = gtk_stack_new();
	options_texttagtable = gtk_text_tag_table_new();
	options_static_texttag = gtk_text_tag_new("static text");
	g_object_set(options_static_texttag, "editable-set", 1, "editable", 0, NULL);

	gtk_text_tag_table_add(options_texttagtable, options_static_texttag);

	drawing_area = gtk_drawing_area_new();
	gtk_widget_set_size_request(drawing_area, 600,500);
	g_signal_connect(G_OBJECT(drawing_area), "draw", G_CALLBACK(draw_callback), user_data);
	gtk_widget_add_events(drawing_area, GDK_BUTTON_PRESS_MASK);
	
	fractal_selector_menu = gtk_menu_new();  
	canopy_menu_item = gtk_menu_item_new_with_label("tree canopy");
	koch_curve_menu_item = gtk_menu_item_new_with_label("koch curve");
	mandelbrot_menu_item = gtk_menu_item_new_with_label("mandelbrot set");
	gtk_menu_attach(fractal_selector_menu,
			canopy_menu_item,
			0, 1, 0, 1);
	gtk_menu_attach(fractal_selector_menu,
			koch_curve_menu_item,
			0, 1, 1, 2);
	gtk_menu_attach(fractal_selector_menu,
			mandelbrot_menu_item,
			0, 1, 2, 3);


	g_signal_connect(koch_curve_menu_item, "activate", G_CALLBACK(koch_curve_menu_item_activate), options_stack);
	g_signal_connect(canopy_menu_item, "activate", G_CALLBACK(canopy_menu_item_activate), options_stack);
	g_signal_connect(mandelbrot_menu_item, "activate", G_CALLBACK(mandelbrot_menu_item_activate), options_stack);
	fractal_selector_menuitem = gtk_menu_item_new_with_label("Select fractal");
	gtk_menu_item_set_submenu(fractal_selector_menuitem, fractal_selector_menu);

	gtk_container_add(GTK_CONTAINER(menu_bar), fractal_selector_menuitem);	
	

	// option grids
	mandelbrot_options_grid = gtk_grid_new();
	iterations_specifier = gtk_entry_new();
	gtk_entry_set_max_length(iterations_specifier, 0);
	gtk_entry_set_placeholder_text(iterations_specifier, "# of iterations");
	significand_length_specifier = gtk_entry_new();
	gtk_entry_set_max_length(significand_length_specifier, 0);
	gtk_entry_set_placeholder_text(significand_length_specifier, "Fraction bit length");
	threshold_specifier = gtk_entry_new();
	gtk_entry_set_max_length(threshold_specifier, 0);
	gtk_entry_set_placeholder_text(threshold_specifier, "threshold");
	mandelbrot_start_button = gtk_button_new_with_label("Draw mandelbrot set");
	mandelbrot_text_buffer = gtk_text_buffer_new(options_texttagtable);
	char *mandelbrot_guidance = "Threshold must be greater\n"
		                    "than 3.0. Once a mandelbrot\n"
	             	            "set image is displayed you\n"
				    "can use the mouse pointer\n"
				    "to select (via clicking on)\n" 
				    "the top left corner and\n"
				    "the bottom right corner of\n"
				    "a rectangle that you wish to\n"
				    "zoom in on. After you have\n"
				    "selected a rectangle to zoom\n"
				    "in on, wait for the machine\n"
				    "to to the calculations and\n"
				    "display the resulting image.\n"
				    "To start with to get a quick\n"
				    "picture you can use iterations\n"
				    "set to 100. fraction bit length\n"
				    "53 (double precision) and\n"
				    "threshold can be set to 50.0";
	gtk_text_buffer_set_text(mandelbrot_text_buffer, mandelbrot_guidance, -1);
	GtkTextIter start_iter;
	GtkTextIter end_iter;
	gtk_text_buffer_get_start_iter(mandelbrot_text_buffer, &start_iter);
	gtk_text_buffer_get_end_iter(mandelbrot_text_buffer, &end_iter);
	gtk_text_buffer_apply_tag_by_name(mandelbrot_text_buffer, "static text", &start_iter, &end_iter);
	mandelbrot_text_view = gtk_text_view_new_with_buffer(mandelbrot_text_buffer);
	g_signal_connect(mandelbrot_start_button, "pressed", G_CALLBACK(check_mandelbrot_options), options_stack);
	g_signal_connect(mandelbrot_start_button, "released", G_CALLBACK(start_mandelbrot_draw), drawing_area);
	gtk_grid_attach(mandelbrot_options_grid, iterations_specifier, 0, 0, 1, 1);
	gtk_grid_attach(mandelbrot_options_grid, significand_length_specifier, 0, 1, 1, 1);
	gtk_grid_attach(mandelbrot_options_grid, threshold_specifier, 0, 2, 1, 1);
	gtk_grid_attach(mandelbrot_options_grid, mandelbrot_start_button, 0, 3, 1, 1);
	gtk_grid_attach(mandelbrot_options_grid, mandelbrot_text_view, 0, 4, 1, 1);
	gtk_stack_add_titled(options_stack, mandelbrot_options_grid, "mandelbrot options", "mandelbrot");

	tree_canopy_options_grid = gtk_grid_new();
	initial_branch_length_specifier = gtk_entry_new();
	gtk_entry_set_max_length(initial_branch_length_specifier, 0);
	gtk_entry_set_placeholder_text(initial_branch_length_specifier, "trunk length");
	GtkWidget *anglegrid = gtk_grid_new();
	GtkWidget *angle_label = gtk_label_new("angle");
	angle_specifier = gtk_scale_button_new(GTK_ICON_SIZE_MENU, M_PI/20.0, M_PI/3.0, M_PI/100.0, NULL);
	gtk_grid_attach(anglegrid, angle_label, 0, 0, 1, 1);
	gtk_grid_attach(anglegrid, angle_specifier, 1, 0, 1, 1);
	ratio = gtk_entry_new();
	gtk_entry_set_max_length(ratio, 0);
	gtk_entry_set_placeholder_text(ratio, "ratio (0.1)");
	depth = gtk_entry_new();
	gtk_entry_set_max_length(depth, 3);
	gtk_entry_set_placeholder_text(depth, "depth");
	canopy_start_button = gtk_button_new_with_label("Draw canopy");
	g_signal_connect(canopy_start_button, "pressed", G_CALLBACK(check_canopy_options), options_stack);
	g_signal_connect(canopy_start_button, "released", G_CALLBACK(start_draw_canopy), drawing_area);
	canopy_text_buffer =  gtk_text_buffer_new(options_texttagtable);
	char *canopy_guidance = "Trunk length is suggested to be\n"
				"set at 100. Ratio should be a\n"
				"decimal number greater than zero\n"
			        "but less than 1. Depth should be\n"
				"greater than 1 and usually less\n"
				"than 20. Angle can be played with";
	gtk_text_buffer_set_text(canopy_text_buffer, canopy_guidance, -1);
	GtkTextIter start;
	GtkTextIter end;
	gtk_text_buffer_get_start_iter(canopy_text_buffer, &start);
	gtk_text_buffer_get_end_iter(canopy_text_buffer, &end);
	gtk_text_buffer_apply_tag_by_name(canopy_text_buffer, "static text", &start, &end);
	canopy_text_view = gtk_text_view_new_with_buffer(canopy_text_buffer);
	gtk_grid_attach(tree_canopy_options_grid, initial_branch_length_specifier, 0, 0, 1, 1);
	gtk_grid_attach(tree_canopy_options_grid, anglegrid, 0, 1, 1, 1);
	gtk_grid_attach(tree_canopy_options_grid, ratio, 0, 2, 1, 1);
	gtk_grid_attach(tree_canopy_options_grid, depth, 0, 3, 1, 1);
	gtk_grid_attach(tree_canopy_options_grid, canopy_start_button, 0, 4, 1, 1); 
	gtk_grid_attach(tree_canopy_options_grid, canopy_text_view, 0, 5, 1, 1);
	gtk_stack_add_titled(options_stack, tree_canopy_options_grid, "canopy options", "canopy");
	
	koch_curve_options_grid = gtk_grid_new();
	koch_iterations_specifier = gtk_entry_new();
	gtk_entry_set_max_length(koch_iterations_specifier, 0);
	gtk_entry_set_placeholder_text(koch_iterations_specifier, "iterations");
	koch_start_button = gtk_button_new_with_label("Draw koch curve");
	g_signal_connect(koch_start_button, "pressed", G_CALLBACK(check_koch_curve_options), options_stack);
	g_signal_connect(koch_start_button, "released", G_CALLBACK(start_koch_curve_draw), drawing_area);
	koch_curve_text_buffer =  gtk_text_buffer_new(options_texttagtable);
	char *koch_curve_guidance = "Keep iterations below 7\n"
	                            "but greated than 0.";
	gtk_text_buffer_set_text(koch_curve_text_buffer, koch_curve_guidance, -1);
	GtkTextIter start1;
	GtkTextIter end1;
	gtk_text_buffer_get_start_iter(koch_curve_text_buffer, &start1);
        gtk_text_buffer_get_end_iter(koch_curve_text_buffer, &end1);
        gtk_text_buffer_apply_tag_by_name(koch_curve_text_buffer, "static text", &start1, &end1);
	koch_curve_text_view = gtk_text_view_new_with_buffer(koch_curve_text_buffer);
	gtk_grid_attach(koch_curve_options_grid, koch_iterations_specifier, 0, 0, 1, 1);
	gtk_grid_attach(koch_curve_options_grid, koch_start_button, 0, 1, 1, 1);
	gtk_grid_attach(koch_curve_options_grid, koch_curve_text_view, 0, 2, 1, 1);
	gtk_stack_add_titled(options_stack, koch_curve_options_grid, "koch curve options", "koch curve");


  	gtk_grid_attach(main_grid, menu_bar, 0, 0, 2, 1);
	gtk_grid_attach(main_grid, drawing_area, 0, 1, 1, 1);
	gtk_grid_attach(main_grid, options_stack, 1, 1, 1, 1);

	gtk_container_add(GTK_CONTAINER(window), main_grid);
  	gtk_widget_show_all (window);
}

int
main (int    argc,
      char **argv)
{
  GtkApplication *app;
  int status;
  fractal_id = -1; // set to tree canopy
  app = gtk_application_new ("org.gtk.example", G_APPLICATION_FLAGS_NONE);
  g_signal_connect (app, "activate", G_CALLBACK (activate), NULL);
  status = g_application_run (G_APPLICATION (app), argc, argv);
  g_object_unref (app);

  return status;
}
