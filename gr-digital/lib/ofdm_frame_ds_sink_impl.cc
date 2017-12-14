/* -*- c++ -*- */
/* 
 * Copyright 2017 Free Software Foundation, Inc.
 * 
 * This file is part of GNU Radio
 * 
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "ofdm_frame_ds_sink_impl.h"
#include <gnuradio/expj.h>
#include <gnuradio/math.h>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <string>


namespace gr {
  namespace digital {

    #define VERBOSE 0
    #define DR_BPSK 2

    inline void
    ofdm_frame_ds_sink_impl::enter_search()
    {
      if(VERBOSE)
	fprintf(stderr, "@ enter_search\n");

      d_state = STATE_SYNC_SEARCH;
    }

    inline void
    ofdm_frame_ds_sink_impl::enter_have_sync()
    {
      if(VERBOSE)
	fprintf(stderr, "@ enter_have_sync\n");

      d_state = STATE_HAVE_SYNC;

      // clear state of demapper
      // d_byte_offset = 0;
      // d_partial_byte = 0;

      d_header = 0;
      d_headerbytelen_cnt = 0;

      // Resetting PLL
      // d_freq = 0.0;
      // d_phase = 0.0;
      // fill(d_dfe.begin(), d_dfe.end(), gr_complex(1.0,0.0));
    }

    inline void
    ofdm_frame_ds_sink_impl::enter_have_header()
    {
      d_state = STATE_HAVE_HEADER;

      // header consists of two 16-bit shorts in network byte order
      // payload length is lower 12 bits
      // whitener offset is upper 4 bits
      d_packetlen = (d_header >> 16) & 0x0fff;
      d_packet_whitener_offset = (d_header >> 28) & 0x000f;
      d_packetlen_cnt = 0;

      if(VERBOSE)
		fprintf(stderr, "@ enter_have_header (payload_len = %d) (offset = %d)\n",
	
		d_packetlen, d_packet_whitener_offset);
    }

 //    char
 //    ofdm_frame_sink_impl::slicer(const gr_complex x)
 //    {
 //      unsigned int table_size = d_sym_value_out.size();
 //      unsigned int min_index = 0;
 //      float min_euclid_dist = std::norm(x - d_sym_position[0]);
 //      float euclid_dist = 0;

 //      for(unsigned int j = 1; j < table_size; j++){
	// euclid_dist = std::norm(x - d_sym_position[j]);
	// if(euclid_dist < min_euclid_dist){
	//   min_euclid_dist = euclid_dist;
	//   min_index = j;
	// }
 //      }
 //      return d_sym_value_out[min_index];
 //    }

    unsigned int ofdm_frame_ds_sink_impl::bpsk_demapper(const gr_complex *in,
						char *out)
    {	
	int BPSK_CONSTEL[2] = {-1,1};
    	float min_dist = 0;
    	unsigned int min_flag_temp = 0;
    	float min_dist_temp = 0;
    	unsigned int bytes_produced = 0;
      unsigned int bytes_offset = 0;
      unsigned int partial_byte = 0;

      int *mod_data = new int[d_occupied_carriers];

    	//decode set1
    	for(unsigned int i=0;i<d_zeros_on_left;i++){
    		for(unsigned int j=0;j<(unsigned int)pow(DR_BPSK,d_downrate-1);j++){
    			min_dist_temp = std::norm(in[i]-d_aliasMap_d0[i][j]);
    			if(min_dist_temp<min_dist||min_dist==0){
    				min_dist=min_dist_temp;
    				min_flag_temp = j;
    			}
    		}
    		mod_data[d_alias_fft_length-d_zeros_on_left+i] = BPSK_CONSTEL[min_flag_temp%DR_BPSK];
    		mod_data[2*d_alias_fft_length-d_zeros_on_left+i] = BPSK_CONSTEL[(min_flag_temp/DR_BPSK)%DR_BPSK];
    		mod_data[3*d_alias_fft_length-d_zeros_on_left+i] = BPSK_CONSTEL[(min_flag_temp/DR_BPSK/DR_BPSK)%DR_BPSK];
    	}min_dist = 0;
    	//decode set2
    	for(unsigned int i=d_zeros_on_left;i<(d_alias_fft_length-d_zeros_on_left);i++){
    		for(unsigned int j=0;j<(unsigned int)pow(DR_BPSK,d_downrate);j++){
    			min_dist_temp = std::norm(in[i]-d_aliasMap_d0[i][j]);
    			if(min_dist_temp<min_dist||min_dist==0){
    				min_dist=min_dist_temp;
    				min_flag_temp = j;
    			}
    		}
    		mod_data[i-d_zeros_on_left] = BPSK_CONSTEL[min_flag_temp%DR_BPSK];
    		mod_data[1*d_alias_fft_length-d_zeros_on_left+i] = BPSK_CONSTEL[(min_flag_temp/DR_BPSK)%DR_BPSK];
    		mod_data[2*d_alias_fft_length-d_zeros_on_left+i] = BPSK_CONSTEL[(min_flag_temp/DR_BPSK/DR_BPSK)%DR_BPSK];
    		mod_data[3*d_alias_fft_length-d_zeros_on_left+i] = BPSK_CONSTEL[(min_flag_temp/DR_BPSK/DR_BPSK/DR_BPSK)%DR_BPSK];
    	}min_dist=0;
    	//decode set3
    	for(unsigned int i=(d_alias_fft_length-d_zeros_on_left);i<d_alias_fft_length;i++){
    		for(unsigned int j=0;j<(unsigned int)pow(DR_BPSK,d_downrate-1);j++){
    			min_dist_temp = std::norm(in[i]-d_aliasMap_d0[i][j]);
    			if(min_dist_temp<min_dist||min_dist==0){
    				min_dist=min_dist_temp;
    				min_flag_temp = j;
    			}
    		}
    		mod_data[i-d_zeros_on_left] = BPSK_CONSTEL[min_flag_temp%DR_BPSK];
    		mod_data[1*d_alias_fft_length-d_zeros_on_left+i] = BPSK_CONSTEL[(min_flag_temp/DR_BPSK)%DR_BPSK];
    		mod_data[2*d_alias_fft_length-d_zeros_on_left+i] = BPSK_CONSTEL[(min_flag_temp/DR_BPSK/DR_BPSK)%DR_BPSK];
    	}

    	//constellation to bytes
    	for(unsigned int i=0;i<d_subcarrier_map.size();i++){
          if(d_derotated_output != NULL){
            d_derotated_output[i]=mod_data[d_subcarrier_map[i]];
          }
          if(bytes_offset==8){
            out[bytes_produced++]=partial_byte;
            bytes_offset=0;
            partial_byte=0;
          }else{
            if(mod_data[d_subcarrier_map[i]]==1){
              partial_byte |= 1 <<(bytes_offset);
            }else{
              partial_byte |= 0 <<(bytes_offset);
            }bytes_offset+=1;
          }
    	}
      delete []mod_data;
    	return bytes_produced;
    }

    void ofdm_frame_ds_sink_impl::generate_bpsk_map(const gr_complex *H_f_in){
	int BPSK_CONSTEL[2] = {-1,1};
    	if(d_downrate==4){
    	//set1
    	for(unsigned int i=0;i<DR_BPSK;i++){
    		for(unsigned int j=0;j<DR_BPSK;j++){
    			for(unsigned int k=0;k<DR_BPSK;k++){
    				//d_idMap_d0[0][i*(unsigned int)pow(DR_BPSK,2)+j*(unsigned int)pow(DR_BPSK,1)+k*(unsigned int)pow(DR_BPSK,0)] = 
    				for(unsigned int l=0;l<d_zeros_on_left;l++){
    					d_aliasMap_d0[l][i*(unsigned int)pow(DR_BPSK,2)+j*(unsigned int)pow(DR_BPSK,1)+k*(unsigned int)pow(DR_BPSK,0)] = H_f_in[d_alias_fft_length+l]*(gr_complex)BPSK_CONSTEL[i]
    					+H_f_in[2*d_alias_fft_length+l]*(gr_complex)BPSK_CONSTEL[j]+H_f_in[3*d_alias_fft_length+l]*(gr_complex)BPSK_CONSTEL[k];
    				}
    			}
    		}
    	}
    	//set2
    	for(unsigned int i=0;i<DR_BPSK;i++){
    		for(unsigned int j=0;j<DR_BPSK;j++){
    			for(unsigned int k=0;k<DR_BPSK;k++){
    				for(unsigned int m=0;m<DR_BPSK;m++){
    				//d_idMap_d0[0][i*(unsigned int)pow(DR_BPSK,2)+j*(unsigned int)pow(DR_BPSK,1)+k*(unsigned int)pow(DR_BPSK,0)] = 
    				for(unsigned int l=d_zeros_on_left;l<(d_alias_fft_length-d_zeros_on_left);l++){
    					d_aliasMap_d0[l][i*(unsigned int)pow(DR_BPSK,3)+j*(unsigned int)pow(DR_BPSK,2)+k*(unsigned int)pow(DR_BPSK,1)+m*(unsigned int)pow(DR_BPSK,0)] = H_f_in[l]*(gr_complex)BPSK_CONSTEL[i]
    					+H_f_in[1*d_alias_fft_length+l]*(gr_complex)BPSK_CONSTEL[j]+H_f_in[2*d_alias_fft_length+l]*(gr_complex)BPSK_CONSTEL[k]
    					+H_f_in[3*d_alias_fft_length+l]*(gr_complex)BPSK_CONSTEL[m];
    				}
    				}
    			}
    		}
    	}
    	//set3
    	for(unsigned int i=0;i<DR_BPSK;i++){
    		for(unsigned int j=0;j<DR_BPSK;j++){
    			for(unsigned int k=0;k<DR_BPSK;k++){
    				//d_idMap_d0[0][i*(unsigned int)pow(DR_BPSK,2)+j*(unsigned int)pow(DR_BPSK,1)+k*(unsigned int)pow(DR_BPSK,0)] = 
    				for(unsigned int l=(d_alias_fft_length-d_zeros_on_left);l<d_alias_fft_length;l++){
    					d_aliasMap_d0[l][i*(unsigned int)pow(DR_BPSK,2)+j*(unsigned int)pow(DR_BPSK,1)+k*(unsigned int)pow(DR_BPSK,0)] = H_f_in[l]*(gr_complex)BPSK_CONSTEL[i]+H_f_in[1*d_alias_fft_length+l]
    					*(gr_complex)BPSK_CONSTEL[j]+H_f_in[2*d_alias_fft_length+l]*(gr_complex)BPSK_CONSTEL[k];
    				}
    			}
    		}
    	}}else{

    	}

    } 

    ofdm_frame_ds_sink::sptr
    ofdm_frame_ds_sink::make(unsigned int fft_length,msg_queue::sptr target_queue,int occupied_carriers,std::string modulation,unsigned int ndelay,unsigned int downrate)
    {
      return gnuradio::get_initial_sptr
        (new ofdm_frame_ds_sink_impl(fft_length, target_queue, occupied_carriers, modulation, ndelay, downrate));
    }

    /*
     * The private constructor
     */
    ofdm_frame_ds_sink_impl::ofdm_frame_ds_sink_impl(unsigned int fft_length,msg_queue::sptr target_queue,int occupied_carriers,std::string modulation,unsigned int ndelay,unsigned int downrate)
      : gr::sync_block("ofdm_frame_ds_sink",
              gr::io_signature::make3(3, 3, sizeof(gr_complex)*fft_length, sizeof(char),sizeof(gr_complex)*fft_length),
              gr::io_signature::make(1, 1, sizeof(gr_complex)*occupied_carriers)),
	d_target_queue(target_queue), d_occupied_carriers(occupied_carriers),
        d_downrate(downrate),d_fft_length(fft_length)
    {
GR_LOG_WARN(d_logger, "The gr::digital::ofdm_frame_sync block has been deprecated.");
      d_alias_fft_length = d_fft_length/d_downrate;
	d_zeros_on_left = ceil((d_fft_length-d_occupied_carriers)/2.0);

      	d_aliasMap_d0 = new gr_complex*[d_alias_fft_length];
	for(unsigned int i=0;i<d_alias_fft_length;i++){
	    d_aliasMap_d0[i]= new gr_complex[(unsigned int)pow(DR_BPSK,d_downrate)];
	}
      	//d_aliasMap_d1 = new gr_complex[d_alias_fft_length];
      	//d_aliasMap_d2 = new gr_complex[d_alias_fft_length];
      	//d_aliasMap_d3 = new gr_complex[d_alias_fft_length];
      	//d_idMap_d0 = new int[3][(unsigned int)pow(DR_BPSK,d_downrate)];


      std::string carriers = "FE7F";

      // A bit hacky to fill out carriers to occupied_carriers length
      int diff = (d_occupied_carriers - 4*carriers.length());
      while(diff > 7) {
		carriers.insert(0, "f");
		carriers.insert(carriers.length(), "f");
		diff -= 8;
      }

      // if there's extras left to be processed
      // divide remaining to put on either side of current map
      // all of this is done to stick with the concept of a carrier map string that
      // can be later passed by the user, even though it'd be cleaner to just do this
      // on the carrier map itself
      int diff_left=0;
      int diff_right=0;

      // dictionary to convert from integers to ascii hex representation
      char abc[16] = {'0', '1', '2', '3', '4', '5', '6', '7',
		      '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};
      if(diff > 0) {
		char c[2] = {0,0};

		diff_left = (int)ceil((float)diff/2.0f);  // number of carriers to put on the left side
		c[0] = abc[(1 << diff_left) - 1];         // convert to bits and move to ASCI integer
		carriers.insert(0, c);

		diff_right = diff - diff_left;	      // number of carriers to put on the right side
		c[0] = abc[0xF^((1 << diff_right) - 1)];  // convert to bits and move to ASCI integer
		carriers.insert(carriers.length(), c);
      }

      // It seemed like such a good idea at the time...
      // because we are only dealing with the occupied_carriers
      // at this point, the diff_left in the following compensates
      // for any offset from the 0th carrier introduced
      int i;
      unsigned int j,k;
      for(i = 0; i < (d_occupied_carriers/4)+diff_left; i++) {
		char c = carriers[i];
		for(j = 0; j < 4; j++) {
	  		k = (strtol(&c, NULL, 16) >> (3-j)) & 0x1;
	  		if(k) {
	    		d_subcarrier_map.push_back(4*i + j - diff_left);
	  		}
		}
      }

      // make sure we stay in the limit currently imposed by the occupied_carriers
      if(d_subcarrier_map.size() > (size_t)d_occupied_carriers) {
		throw std::invalid_argument("ofdm_frame_sink_impl: subcarriers allocated exceeds size of occupied carriers");
      }

      d_bytes_out = new char[d_occupied_carriers];
      //d_dfe.resize(occupied_carriers);
    //fill(d_dfe.begin(), d_dfe.end(), gr_complex(1.0,0.0));

      //set_sym_value_out(sym_position, sym_value_out);

      enter_search();
     }

    /*
     * Our virtual destructor.
     */
    ofdm_frame_ds_sink_impl::~ofdm_frame_ds_sink_impl()
    {
	delete [] d_bytes_out;
	for(unsigned int i=0;i<d_alias_fft_length;i++){
	   delete []d_aliasMap_d0[i];	
	}
	delete []d_aliasMap_d0;
    //  delete [] d_aliasMap_d1;
    //  delete [] d_aliasMap_d2;
    //  delete [] d_aliasMap_d3;
    }

    int
    ofdm_frame_ds_sink_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex*)input_items[0];
      const char *sig = (const char*)input_items[1];
      const gr_complex *H_f_in = (const gr_complex*)input_items[2];
      // If the output is connected, send it the derotated symbols
      if(output_items.size() >= 1)
		    d_derotated_output = (gr_complex *)output_items[0];
      else
		    d_derotated_output = NULL;

      	gr_complex *in_d0 = new gr_complex[noutput_items*d_alias_fft_length];
      	// gr_complex *in_d1 = new gr_complex[noutput_items*d_alias_fft_length];
      	// gr_complex *in_d2 = new gr_complex[noutput_items*d_alias_fft_length];
      	// gr_complex *in_d3 = new gr_complex[noutput_items*d_alias_fft_length];


      for(int i=0;i<noutput_items;i++){
      	for(unsigned int j=0;j<d_alias_fft_length;j++){
      		for(unsigned int k=0;k<d_downrate;k++){
      			in_d0[i*d_fft_length+j] +=in[i*d_fft_length+k*d_alias_fft_length+j]
      			*gr_expj(2*M_PI*(k*d_alias_fft_length+j+1)*DELAY0/d_fft_length);
      			// in_d1[i*d_fft_length+j] +=in[i*d_fft_length+k*d_alias_fft_length+j]
      			// *gr_expj(2*M_PI*(k*d_alias_fft_length+j+1)*DELAY1/d_fft_length);
      			// in_d2[i*d_fft_length+j] +=in[i*d_fft_length+k*d_alias_fft_length+j]
      			// *gr_expj(2*M_PI*(k*d_alias_fft_length+j+1)*DELAY2/d_fft_length);
      			// in_d3[i*d_fft_length+j] +=in[i*d_fft_length+k*d_alias_fft_length+j]
      			// *gr_expj(2*M_PI*(k*d_alias_fft_length+j+1)*DELAY3/d_fft_length);
      		}
      	}
      }

      unsigned int j = 0;
      unsigned int bytes=0;

      switch(d_state){
      	case STATE_SYNC_SEARCH:
      		if(sig[0]){
            if(d_downrate!=0)
              generate_bpsk_map(&H_f_in[0]);
      			enter_have_sync();
      		}
      		break;

      	case STATE_HAVE_SYNC:
      		bytes = bpsk_demapper(&in_d0[0],d_bytes_out);

      		j = 0;
			    while(j < bytes) {
	  			  d_header = (d_header << 8) | (d_bytes_out[j] & 0xFF);
	  				j++;

	  			  if(++d_headerbytelen_cnt == HEADERBYTELEN) {
	    		   // we have a full header, check to see if it has been received properly
	    			if(header_ok()) {
	      				enter_have_header();
	      			 	while((j < bytes) && (d_packetlen_cnt < d_packetlen)) {
	      			 		d_packet[d_packetlen_cnt++] = d_bytes_out[j++];
	      			 	}

	      				if(d_packetlen_cnt == d_packetlen) {
							message::sptr msg =message::make(0, d_packet_whitener_offset, 0, d_packetlen);
							memcpy(msg->msg(), d_packet, d_packetlen_cnt);
							d_target_queue->insert_tail(msg);		// send it
							msg.reset();  				// free it up

							enter_search();
						}
					}
					else {
						enter_search();				// bad header
					}
				}	
			}
			break;


      	case STATE_HAVE_HEADER:
      		bytes = bpsk_demapper(&in_d0[0],d_bytes_out);

      		j = 0;
			while(j < bytes) {
	  			d_packet[d_packetlen_cnt++] = d_bytes_out[j++];

	  			if (d_packetlen_cnt == d_packetlen){		// packet is filled
	    		// build a message
	    		// NOTE: passing header field as arg1 is not scalable
	    			message::sptr msg =message::make(0, d_packet_whitener_offset, 0, d_packetlen_cnt);
	   				memcpy(msg->msg(), d_packet, d_packetlen_cnt);

	   				d_target_queue->insert_tail(msg);		// send it
	    			msg.reset();  				// free it up

	    			enter_search();
	    			break;
	  			}
			}
			break;

		default:
			assert(0);
		} // switch

    delete []in_d0;
		return 1;
    }

  } /* namespace digital */
} /* namespace gr */

