cd /tmp && wget https://www.python.org/ftp/python/2.7/Python-2.7.tgz && tar -xzf Python-2.7.tgz && cd Python-2.7 && CC=/usr/bin/clang ./configure && make &&  sudo make altinstall
cd /tmp && git clone https://github.com/numpy/numpy && cd numpy && git checkout v1.13.3 && git submodule update --init --recursive  && sudo python setup.py install
