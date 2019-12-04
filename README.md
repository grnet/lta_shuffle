# lta_shuffle
Repository for the implementation of A Non-Interactive Shuffle Argument With Low Trust Assumptions

Pre-requisites
--------------
On Ubuntu/Debian do:

```
apt install python3 python3-dev libffi-dev libssl-dev build-essential
```

Python 3.5+ required.


Install requirements
--------------------

```
pip install -r requirements.txt
```


Multi-Party Computation Demo
----------------------------
To run a multi-party computation use the following for each peer:

```
python demo.py [-h] [--fixed] [--nb] --name name [--config config] [N]
```

where both peers should use the same value of N > 1 (see below).

For example, using the sample config with two peers, run

```
python demo.py --name one --config bbparty.config.local
```

```
python demo.py --name two --config bbparty.config.local
```

from different terminals.


Secure values of N
------------------
In the current implementation (based upon `bplib`), in order to produce secure

results compliant with the paper, N should be such that

`N + 1 | p - 1`

where

`p = 16798108731015832284940804142231733909759579603404752749028378864165570215949`

is the prime used by OpenPairing (https://github.com/dfaranha/OpenPairing).

Some such values are

`2`, `3`, `5`, `6`, `11`, `13`, `20`, `27`, `41`, `83`, `640`, `1281`, `1922`, `2563`, `3845`, `4486`, `7691`, ...

NOTE: Avoid `N > 27`, since local computations `(O(N ^ 2))` are not yet optimized
and the execution becomes very slow for large values of `N`.


Contributors
------------

Giorgos Korfiatis

Antonis Angelakis

Giorgos Tsoukalas

Foteinos Mergoupis
