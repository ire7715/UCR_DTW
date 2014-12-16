<h1>UCR_DTW</h1>

&emsp;An implementation of "<b>Searching and Mining Trillions of Time Series Subsequences under Dynamic Time Warping, 2011</b>"
 
-------

<h2>Improvement</h2>
 - Fix the sorting bug. (Origin version wouldn't swap when the difference of  two values is less than 1.)
 - Delay the z-normalize for t sequence. (Only dtw apply the tz array, so there is no need to do z-normalize when keogh_data)
 - <b>Jump pruning</b>: When we are pruning by sorted lower bound, we can find out that which index is the most front index(j). And instead of compute the next linearly, we jump to (j+1) directly. Jump pruning prunes over 50% for query2.txt.

-------

<h4>Result without jump pruning</h4>
```
>ruby UCR_DTW.rb Data.txt Query2.txt 0.1 -nj
Location: 430264
Distance: 3.7906991699901917
Data scanned: 1000000
Execution time: 121.967059865

Pruned by Jump: 0.0%
Pruned by LB_Kim: 22.6827%
Pruned by LB_Keogh: 35.9776%
Pruned by LB_Keogh2: 36.5951%
DTW Calcuation  : 4.744599999999991
```

<h4>Result with jump pruning</h4>
```
>ruby UCR_DTW.rb Data.txt Query2.txt 0.1
Location: 430264
Distance: 3.7906991699901917
Data scanned: 1000000
Execution time: 59.746849525
 
Pruned by Jump: 57.7814%
Pruned by LB_Kim: 16.4837%
Pruned by LB_Keogh: 15.1936%
Pruned by LB_Keogh2: 7.7391%
DTW Calcuation  : 2.802199999999999
```
