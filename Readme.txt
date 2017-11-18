# Author:HZL Time:2017/11/17 Location:Room W316, South Building 1, Huazhong University of Science and Technology
 
This gnuradio project is used to reduce the power in WiFi communication,and it is based on the paper:"Sampleless Wi-Fi: Bringing Low Power to Wi-Fi Communications.Wei Wang, Yingjie Chen, Lu Wang, Qian Zhang.IEEE/ACM Trans. Netw., vol.25, no.3, 2017."

Here are some tips for correctly using this project:

1.ofdm_oversampler block in gr-digital file and demod_vcc block in gr-fft block are totally new c++ block which belong to digital module and fft module.You should use gr_modtool to add these blocks into digital and fft module(eg. $ gr_modtool add ofdm_oversampler).What's more,you can get the argument list which needed for creating a new block in header files(eg. ofdm_sampler.h).

When you successfully created these blocks,you can directly copy these files in our project into your new files.

2.Other files are based on the original files in gnuradio,so you directly copy these files in our project into your gnuradio files.

3.Finally,you should rebuild the gnuradio using the commands below:
$ cd ../gnuradio-3.7.11/build (go to the build file of gnuradio)
$ cmake ../
$make
$sudo make install
After finishing these steps,you can test it in your USRP platform.


