from gcc:9.4

ADD SAINTexpress_v3.6.3__2018-03-09.tar.gz saint.source

RUN sed "s/setprecision(.)/setprecision(8)/" /saint.source/SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-int/main.hpp > /saint.source/SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-int/newmain.hpp; mv /saint.source/SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-int/newmain.hpp /saint.source/SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-int/main.hpp
RUN sed "s/setprecision(.)/setprecision(8)/" /saint.source/SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-spc/main.hpp > /saint.source/SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-spc/newmain.hpp; mv /saint.source/SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-spc/newmain.hpp /saint.source/SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-spc/main.hpp


RUN make -C /saint.source/SAINTexpress_v3.6.3__2018-03-09/ -j

RUN ln -s /saint.source/SAINTexpress_v3.6.3__2018-03-09/bin/SAINTexpress-spc  /bin/
RUN ln -s /saint.source/SAINTexpress_v3.6.3__2018-03-09/bin/SAINTexpress-int  /bin/
