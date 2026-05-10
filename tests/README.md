To test the different functionalities of `tinyMapper.sh`:

```sh
# ChIP
./tinyMapper.sh -m ChIP -s tests/testChIP.IP -g tests/R64-1-1/R64-1-1
./tinyMapper.sh -m ChIP -s tests/testChIP.IP -i tests/testChIP.input -g tests/R64-1-1/R64-1-1
./tinyMapper.sh -m ChIP -s tests/testChIP.IP -i tests/testChIP.input -g tests/R64-1-1/R64-1-1 -c tests/CBS138/CBS138

# RNA
./tinyMapper.sh -m RNA -s tests/testRNA -g tests/R64-1-1/R64-1-1

# MNase
./tinyMapper.sh -m MNase -s tests/testMNase -g tests/R64-1-1/R64-1-1

# HiC
./tinyMapper.sh -m HiC -s tests/testHiC -g tests/R64-1-1/R64-1-1 --resolutions 1000,2000,8000 --restriction 'DpnII,HinfI'

# ATAC
## tbd...
```
