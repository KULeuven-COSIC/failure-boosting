# failure-boosting

This code allows to calculate the cost of (directional) failure boosting attacks as discussed in https://eprint.iacr.org/2021/193.pdf

* main_multitarget_limited.py: calculates the cost of a multitarget attack
* main_attackcostvstargets.py: calculates the cost of multitarget attack in function of the number of targets

Generation of failure boosting curves is costly, but should be done only once (intermediate results are saved). As such running the main functions for the first time will take approximately 1 day to 1 week of time per scheme and per set of constraints. If you would like the intermediate results you can contact the authors.
