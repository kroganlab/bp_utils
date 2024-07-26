from gcc:9.4

ADD SAINTexpress_v3.6.3__2018-03-09.tar.gz saint.source

RUN make -C /saint.source/SAINTexpress_v3.6.3__2018-03-09/ -j

RUN ln -s /saint.source/SAINTexpress_v3.6.3__2018-03-09/bin/SAINTexpress-spc  /bin/
RUN ln -s /saint.source/SAINTexpress_v3.6.3__2018-03-09/bin/SAINTexpress-int  /bin/
