# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import sys
import numpy as np
from datetime import date, datetime, timedelta
from casacore.quanta import quantity
from casacore.measures import measures


# %%
def rise_set_UTC(sourcename, day_UT_py, obspos):
    """
    Calculate rise an set UTC times of a source on a day

    Examples
    --------
    >>> from datetime import date
    >>> from casacore.measures import measures
    >>> rs_utc_py=rise_set_UTC('SUN', date(2020,11,28),
    ...                        measures().observatory('WSRT'))
    {'rise': datetime.datetime(2020, 11, 28, 8, 9, 6, 42707),
     'set': datetime.datetime(2020, 11, 28, 14, 30, 10, 217090)}
    """
    dm=measures()
    # Set frame:
    #   set position via ITRF position
    obspos_ITRF = dm.position('ITRF', v0=str(obspos[0])+'m',
                              v1=str(obspos[1])+'m', v2=str(obspos[2])+'m')
    dm.doframe(obspos_ITRF)
    #   set epoch via python datetime
    if type(day_UT_py) is datetime:
        # Remove time (keep only date)
        day_UT_py = day_UT_py.date()
    day_UT = quantity(day_UT_py.isoformat())
    day_UT_epoch = dm.epoch('UTC', day_UT)
    dm.doframe(day_UT_epoch)
    # Let casacore calculate rise and set times
    #   (rise() method returns both)
    #   which is given in LST as hour angle in deg.
    rise_set = dm.rise(dm.direction(sourcename), ev='5deg')
    # Determine days
    day_lst_epoch = dm.measure(day_UT_epoch, 'LMST')
    day_lst_qnt = quantity(day_lst_epoch['m0'])
    day_lst_qnt.set_value(int(day_lst_qnt.get_value()))
    rs_utc_py = {}
    # Loop over rise and set
    for rsi in rise_set:
        rs_lst_qnt = rise_set[rsi].to_time() + day_lst_qnt
        rs_lst_epoch = dm.epoch('LMST', rs_lst_qnt)
        rs_utc_epoch = dm.measure(rs_lst_epoch, 'UTC')
        rs_utc_py[rsi] = datetime.utcfromtimestamp(quantity(rs_utc_epoch['m0'])\
            .to_unix_time())
    return rs_utc_py


# %%
from ilisa.antennameta.antennafieldlib import getArrayBandParams
stnpos, _rot, _relpos, _intilepos = getArrayBandParams('UK902','HBA')

# %%
margin = timedelta(minutes=60)  # time margin at end of a scan
one_day = timedelta(days=1)
if len(sys.argv)>1:
    ref_day0 = date.today()
else:
    ref_day0 = date.today()+one_day
#tomorrow=date.today()+one_day
day1 = ref_day0+one_day
stnpos_ar=np.array(stnpos).squeeze()
rs_0 = rise_set_UTC('SUN', ref_day0, stnpos_ar)
rs_1 = rise_set_UTC('SUN', day1, stnpos_ar)
print("Sunrise start: {:%Y-%m-%dT%H:%M:%S},".format(rs_0['rise'])
      +" duration: {:.0f}s, againafter: {:.0f}+{}".format(
       (rs_0['set']-rs_0['rise']).total_seconds(),
       (rs_1['rise']-rs_0['set']-margin).total_seconds(),
       margin.total_seconds()))
# %%
