Readme
================
Christopher Dory
2025-12-16

# Purpose

The purpose of this repository is to hold two functions,
“calculate_stream_depletions” and “map_stream_depletions”. These
functions estimate by analytical equations the depletions experienced on
stream reaches as a function of pumping and a map of the impact of wells
on nearby streams respectively. The former uses linear superposition,
assuming the aquifer responds linearly, to overlay the pumping at well
locations and can accomodate user specified variable pumping rates. The
latter assumes a unit constant pumping rate to estimate at a given well
location how much of the extracted water is from storage and how much
from depletion after user specified timesteps.The author would like to
stress that there are no new ideas contained within this code, and they
serve rather as convienent ways for the user to create depletion maps
and estimates without writing the methods themselves. Please see within
each folder more detailed documentation and references.
