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

#ifndef INCLUDED_DIGITAL_OFDM_FRAME_DS_SINK_IMPL_H
#define INCLUDED_DIGITAL_OFDM_FRAME_DS_SINK_IMPL_H

#include <gnuradio/digital/ofdm_frame_ds_sink.h>

namespace gr {
  namespace digital {

    class ofdm_frame_ds_sink_impl : public ofdm_frame_ds_sink
    {
     private:
      enum state_t {STATE_SYNC_SEARCH, STATE_HAVE_SYNC, STATE_HAVE_HEADER};

      static const int MAX_PKT_LEN    = 4096;
      static const int HEADERBYTELEN   = 4;
      static const int DELAY0 = 0;
      int BPSK_CONSTEL[2] = {-1,1};

      msg_queue::sptr  d_target_queue;	// where to send the packet when received
      state_t            d_state;
      unsigned int       d_header;		// header bits
      int                d_headerbytelen_cnt;	// how many so far

      char              *d_bytes_out;           // hold the current bytes produced by the demapper    

      int                d_occupied_carriers;
      unsigned int       d_alias_fft_length;
      unsigned int       d_downrate;
      unsigned int       d_fft_length;
      unsigned int       d_zeros_on_left;

      char               d_packet[MAX_PKT_LEN];    // assembled payload
      int 	         d_packetlen;              // length of packet
      int                d_packet_whitener_offset; // offset into whitener string to use
      int		 d_packetlen_cnt;          // how many so far

      gr_complex * d_derotated_output;  // Pointer to output stream to send deroated symbols out
      gr_complex ** d_aliasMap_d0;

      // std::vector<gr_complex>    d_sym_position;
      // std::vector<char>          d_sym_value_out;
      // std::vector<gr_complex>    d_dfe;
      // unsigned int d_nbits;

      // char d_resid;
      // unsigned int d_nresid;
      // float d_phase;
      // float d_freq;
      // float d_phase_gain;
      // float d_freq_gain;
      // float d_eq_gain;

      std::vector<int> d_subcarrier_map;

     protected:
      void enter_search();
      void enter_have_sync();
      void enter_have_header();
  
      bool header_ok()
      {
	// confirm that two copies of header info are identical
	return ((d_header >> 16) ^ (d_header & 0xffff)) == 0;
      }
  
      // char slicer(const gr_complex x);

      unsigned int bpsk_demapper(const gr_complex *in,
            char *out);
      void generate_bpsk_map(const gr_complex *H_f_in);

      // bool set_sym_value_out(const std::vector<gr_complex> &sym_position, 
			   //   const std::vector<char> &sym_value_out);

     public:
      ofdm_frame_ds_sink_impl(unsigned int fft_length,msg_queue::sptr target_queue,int occupied_carriers,std::string modulation,unsigned int ndelay,unsigned int downrate);
      ~ofdm_frame_ds_sink_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_OFDM_FRAME_DS_SINK_IMPL_H */

