# Rank-based Gene Pair Analysis (rGPA)
Rank-based gene pair analysis

## Usage
how to use.

## Install 
```
module XOR(output out1,  input in1, in2);
  wire w1, w2, w3, w4;
  not (w1, in1);
  not (w2, in2);
  not (w3, in1, w2);
  not (w4, in2, w1);
  or (out1, w3, w4);
endmodule
```


## Contact Us
If you would like to receive updates from the Cello team regarding bug fixes, patches, feature updates or if you would like to contact the Cello team, please check the links in [CONTACT.md](CONTACT.md)


