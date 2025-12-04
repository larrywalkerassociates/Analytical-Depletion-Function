Arguments
================
Christopher Dory
2025-11-18

# Variable documentation

<b><font size = "3">streams</font></b>: *sf object*, series of
linestrings describing the streams in the domain. <br/> <br/>
<b><font size = "3">streams_are_points</font></b>: *boolean* default
*FALSE*, describes whether the streams being passed are linestrings, in
which the program will transform them to points for calculations, or
whether they are pre-processed by the user to points. An argument of
*FALSE* equates to telling the program you are passing linestrings.
<br/> <br/> <b><font size = "3">stream_id_key</font></b>: *char* default
*NULL*, if streams are already points this argument describes an ID
column in that point set. The ID describes what reach each point belongs
to. IDs must be sequentially numeric (i.e. 1,2,3,4…). <br/> <br/>
<b><font size = "3">wells</font></b>: *sf object*, series of points that
represent the wells in the domain. <br/> <br/>
<b><font size = "3">wells_id_key</font></b>: *char* default *NULL*,
describes a column in that point set where an ID can be found in the
wells data. The ID describes what reach each point belongs to. IDs must
be sequentially numeric (i.e. 1,2,3,4…). <br/> <br/>
<b><font size = "3">pumping</font></b>: *matrix*, a matrix where the
number of rows equals the number of wells, and the number of columns
equals the number of timesteps. Each row, column pair contains a pumping
rate for that timestep. For example, if the timesteps are in days and
m\[1,1\] = 100 this is equivalent to a pumping of 100 units on day 1 for
well 1. <br/> <br/> <b><font size = "3">subwatersheds</font></b>: *sf
object* default *NULL*, if the proximity criteria is set to ‘adjacent’
or ‘adjacent+expanding’ this argument is necessary. It describes what
subwatersheds are to be used in the proximity criteria assignment. For
each well streams within the subwatershed of the well or adjacent
subwatersheds are considered to be effected by the wells pumping. <br/>
<br/> <b><font size = "3">influence_radius</font></b>: *numeric* default
*NULL*, if the proximity criteria is set to ‘local area’, ‘expanding’,
or ‘adjacent+expanding’ this argument is necessary. It describes the
radius around each well that they influence. Streams within this radius
are considered to be effected by the wells pumping. <br/> <br/>
<b><font size = "3">proximity_criteria</font></b>: *char* default
*‘local area’*, describes what proximity criteria to use when assigning
what streams are effected by the pumping of a given well. Can take
arguments of *‘local area’, ‘whole domain’, ‘adjacent’,
‘adjacent+expanding’, ‘expanding’*. In this program *expanding* and
*local area* are equivalent, however the number passed to
*influence_radius* by the user should necessarily be different given
what they represent. <br/> <br/>
<b><font size = "3">apportionment_criteria</font></b>: *char* default
*‘inverse distance’*, describes how to assign for each well what
fraction of its pumping goes to each stream it effects. The sum of the
fractions, no matter the argument passed here, will always total to 1.
Accepts *‘inverse distance’, ‘inverse distance squared’, ‘web’, ‘web
squared’, ‘thiessen polygon’*. <br/> <br/>
<b><font size = "3">geologic_apportionment</font></b>: *boolean* default
*FALSE*, describes whether to scale apportionment by the ratio between
the transmissivity of the geologic unit of the reach and the geologic
unit of the well $\frac{T_{r}}{T_{w}}$. The logic behind this is that a
reach in more transmissive units will receive more depletion than those
in less transmissive units. For example if two reaches are equidistant
from a well, and the reaches are in two different geologic units, they
should not receive equal apportionment. <br/> <br/>
<b><font size = "3">analytical_model</font></b>: *char* default
*‘glover’*, describes what analytical model to use for the calculation
of stream depletions. Accepts *‘glover’, ‘hunt’, ‘hantush’*. *‘Hunt’*
and *‘hantush’* models are equivalent under the special case that
$\lambda = 2\frac{T}{L}$ according to Reeves (2008). <br/> <br/>
<b><font size = "3">depletion_potential_criteria</font></b>: *char*
default *‘fractional’*, describes how to average fractional depletions.
Accepts *‘global’, ‘fractional’,
‘fractional+pumping’*,*‘distance’*,*‘fractional+distance’*,*‘all’*.
*‘Global’* is equivalent to saying at each timestep take the average of
depletion potential for all wells effecting reach ‘r’. *‘Fractional’*
takes a weighted average by frac/mean(frac) from wells ‘w’ to reach ‘r’.
*‘Fractional+pumping’* takes a weighted average by (pumping \*
frac)/mean(pumping \* frac) from wells ‘w’ to reach ‘r’. etc.. <br/>
<br/> <b><font size = "3">sdf_averaging_criteria</font></b>: *char*
default *‘fractional’*, describes what how to average calculated sdf
(units of time). Accepts *‘global’, ‘fractional’, ‘fractional+pumping’*.
*‘Global’* is equivalent to saying at each timestep take the average of
depletion potential for all wells effecting reach ‘r’. *‘Fractional’*
takes a weighted average by frac/mean(frac) from wells ‘w’ to reach ‘r’.
*‘Fractional+pumping’* takes a weighted average by (pumping \*
frac)/mean(pumping \* frac) from wells ‘w’ to reach ‘r’. sdf defined by
Jenkins (1968) as $\frac{d^2S}{T}$. <br/> <br/>
<b><font size = "3">custom_sdf_time</font></b>: *char* default *NULL*,
describes whether to find the time when the depletions in reach ‘r’
equal the pumping in well ‘w’ multiplied by *custom_sdf_time*. Accepts
any number, but recommended to chose between 0.1 and 0.9. If set to
*NULL* no calculations are done. If time to depletions starts between
timestep 0 and 1, -Inf returned. If time to depletions is after final
timestep Inf returned. <br/> <br/>
<b><font size = "3">data_out_dir</font></b>: *char* default *getwd()*,
describes where to place the output data. <br/> <br/>
<b><font size = "3">diag_out_dir</font></b>: *char* default *getwd()*,
describes where to place the diagnostic data (currently limited to
log.txt). <br/> <br/>
<b><font size = "3">suppress_loading_bar</font></b>: *boolean* default
*TRUE*, describes whether to suppress loading bars in console. Value of
*TRUE* does not give loading bars. <br/> <br/>
<b><font size = "3">suppress_console_messages</font></b>: *boolean*
default *TRUE*, describes whether to supress user messages to console.
Value of *TRUE* does not give console messages. <br/> <br/>
<b><font size = "3">stor_coef_key</font></b>: *char* default *‘Stor’*,
describes the name of the column in the well data that stores the
storage coefficient data. <br/> <br/>
<b><font size = "3">well_transmissivity_key</font></b>: *char* default
*‘Tr’*, describes the name of the column in the well data that stores
transmissivity data ($Kb$). <br/> <br/>
<b><font size = "3">stream_transmissivity_key</font></b>: *char* default
*NULL*, describes the name of the column in the stream data that stores
transmissivity data ($Kb$). Only necessary if *geologic_apportionment*
is set to *TRUE*. <br/> <br/>
<b><font size = "3">leakance_key</font></b>: *char* default *NULL*,
describes the name of the column in the stream data that stores leakance
data, where leakance is a product of the conductivity of the aquifer
multiplied by the streambed clogging layer thickness divided by the
streambed clogging layer conductivity $\frac{Kb^\text{'}}{K^\text{'}}$.
Argument is required if *analytical_model* is set to *‘hantush’*.
Definition taken from Reeves (2008) referencing Hantush (1965). <br/>
<br/> <b><font size = "3">lambda_key</font></b>: *char* default *NULL*,
describes the name of the column in the stream data that stores lambda
data, where lambda is a product of the width of the river multiplied by
the conductivity of the streambed clogging layer divided by the
thickness of the streambed clogging layer $w_{r}*\frac{K_{r}}{b_{r}}$.
Argument is required if *analytical_model* is set to *‘hunt’*.
Definition taken from Zipper et al. (2019) referencing Hunt (1999).
<br/> <br/> <b><font size = "3">prec</font></b>: *numeric* default *80*,
describes how many bits of precision to use in the Rmpfr package when
calculating exponential values. Required for the *‘hunt’* and
*‘hantush’* analytical models. Exponential values in the *‘hunt’* and
*‘hantush’* analytical models can often return **Inf** using base
numerics, use of *prec* in Rmpfr alleviates this. Higher *prec* values
will slow down processing. Default taken from streamDepletr package by
Zipper (2019). <br/> <br/> <br/> <br/>

# References

Glover and Balmer (1954) <https://doi.org/10.1029/TR035i003p00468> <br/>
Hantush (1965) <https://doi.org/10.1029/JZ070i012p02829> <br/> Hunt
(1999) <https://doi.org/10.1111/j.1745-6584.1999.tb00962.x> <br/>
Jenkins (1968) <https://pubs.usgs.gov/twri/twri4d1/pdf/twri_4-D1_a.pdf>
<br/> Reeves (2008) <https://doi.org/10.3133/ofr20081166> <br/> Zipper
(2019)
<https://cran.r-project.org/web/packages/streamDepletr/index.html> <br/>
Zipper et al. (2019) <https://doi.org/10.1029/2018WR024403> <br/>
