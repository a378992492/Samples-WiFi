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


#ifndef INCLUDED_DIGITAL_OFDM_FRAME_DS_SINK_H
#define INCLUDED_DIGITAL_OFDM_FRAME_DS_SINK_H

#include <gnuradio/digital/api.h>
#include <gnuradio/sync_block.h>
#include <gnuradio/msg_queue.h>

namespace gr {
  namespace digital {

    /*!
     * \brief <+description of block+>
     * \ingroup digital
     *
     */
    class DIGITAL_API ofdm_frame_ds_sink : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<ofdm_frame_ds_sink> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of digital::ofdm_frame_ds_sink.
       *
       * To avoid accidental use of raw pointers, digital::ofdm_frame_ds_sink's
       * constructor is in a private implementation
       * class. digital::ofdm_frame_ds_sink::make is the public interface for
       * creating new instances.
       */
      static sptr make(unsigned int fft_length,msg_queue::sptr target_queue,int occupied_carriers,std::string modulation,unsigned int ndelay=0,unsigned int downrate=4);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_OFDM_FRAME_DS_SINK_H */

