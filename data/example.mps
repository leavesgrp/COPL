NAME           example
*你猜猜会发生什么事情
ROWS
 L  r1
 L  r2
 L  r3
 N  cost
COLUMNS
    x01       r1                 9.0   r2                 4.0
    x01       r3                 3.0   cost              -7.0
    x02       r1                 4.0   r2                 5.0
    x02       r3                10.0   cost             -12.0
RHS
    rhs       r1               360.0   r2               200.0
    rhs       r3               300.0
BOUNDS
 LO           x01                  1
 UP           x01                100
 LO           x02                  1
 UP           x02                100
ENDATA
