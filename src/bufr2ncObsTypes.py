#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import sys
import re
import netCDF4
from netCDF4 import Dataset
import struct
import datetime as dt

import bufr2ncCommon as cm

############################################################################
# SUBROUTINES
############################################################################

def SplitMsgDate(yyyymmddhh):
    # This routine will take an integer date yyyymmddhh and return the
    # datetime equivalent.
    DateString = str(yyyymmddhh)
    Dtime = dt.datetime(int(DateString[0:4]), int(DateString[4:6]),
                        int(DateString[6:8]), int(DateString[8:10]))

    return Dtime

def MakeDate(Dtime):
    # This routine will take in a datetime object and return an integer date yyyymmdd.
    DateString = "%0.4i"%(Dtime.year) + "%0.2i"%(Dtime.month) + "%0.2i"%(Dtime.day)

    return int(DateString)

def MakeTime(Dtime):
    # This routine will take in a datetime object and return an integer time hhmmss.
    TimeString = "%0.2i"%(Dtime.hour) + "%0.2i"%(Dtime.minute) + "%0.2i"%(Dtime.second)

    return int(TimeString)

def MakeEpochTime(Dtime):
    # This routine will take in a datetime object and return time since the
    # epoch (Jan 1, 1970, 00Z).
    #
    # Need to mark the Dtime object as UTC time before accessing the epoch timestamp.
    EpochTime = Dtime.replace(tzinfo=dt.timezone.utc).timestamp()

    return EpochTime

def DateOffsetToAbsTime(MsgDate, TimeOffset):
    # This routine will calculate the absolute date and time values from
    # input values holding the message date and a time offset.
    #
    # The idea is to convert the message date (integer in the form YYMMDDHH) to a
    # DateTime object, then use timedelta to adjust for the offset.
    #
    # Return the time in two forms
    #    AbsDate: integer in the form of YYYYMMDD
    #    AbsTime: integer in the form of HHMMSS
    #
    #    EpochTime: float, number of seconds from the epoch (Jan 1, 1970, 00Z)

    MsgDtime = SplitMsgDate(MsgDate)

    # At this point it is possible that OdateOffset is a scalar value
    # (indicated by size == 1). If so, then cast OdateOffest into a vector
    # of size one for the following loop.
    # The datetime routines don't accept arrays as arguments so
    # need to calculate absolute times one element at a time.
    Nlevs = TimeOffset.size
    if (Nlevs == 1):
        OdateOffset = np.array([ TimeOffset ])
    else:
        OdateOffset = TimeOffset
    AbsDate = np.empty(Nlevs, dtype=np.int)
    AbsTime = np.empty(Nlevs, dtype=np.int)
    EpochTime = np.empty(Nlevs, dtype=np.float)
    for i in range(Nlevs):
        OffsetDtime = dt.timedelta(hours=float(OdateOffset[i]))
        AbsDtime = MsgDtime + OffsetDtime

        AbsDate[i] = MakeDate(AbsDtime)
        AbsTime[i] = MakeTime(AbsDtime)
        EpochTime[i] = MakeEpochTime(AbsDtime)

    return [ AbsDate, AbsTime, EpochTime ]

def BufrFloatToString(Fval):
    # This routine will interperet a floating point (double) as eight character
    # bytes, thus converting the float value to a string. This is necessary since
    # libbufr routines return string values as doubles.

    # Need to place the value into a numpy array first before converting.
    FvalArray = np.array([ Fval ])

    # Use the struct unpack command to interpret the value as a series of
    # character bytes.
    ByteList = struct.unpack('8c', FvalArray)

    # Join the bytes together and interpret as an ascii string
    Sval = bytes.join(b'', ByteList).decode('ascii').strip()

    return(Sval)

def PrepBufrDictCreateVar(Vname, ObsLat, ObsLon, AbsDate, AbsTime,
                          EpochTime, P, Var, VarErr, VarQc):
    # This routine is intended to be used to create a variable entry in the dictionary
    # of a conventional obs type (radiosonde, aircraft, etc.) when reading
    # from a prepBUFR file.
    VerrName = Vname + "_err"
    VqcName = Vname + "_qc"

    ObsDict = {
       cm.NC_LAT_NAME : [ ],
       cm.NC_LON_NAME : [ ],
       cm.NC_DATE_STAMP_NAME : [ ],
       cm.NC_TIME_STAMP_NAME : [ ],
       cm.NC_TIME_NAME : [ ],
       cm.NC_P_NAME : [ ],
       Vname : [ ],
       VerrName : [ ],
       VqcName : [ ]
       }

    # Walk through arrays and copy values as you go into the dictionary. If there
    # are any missing values, don't do the copy. Missing value from BUFR is 1e11.
    # Call a value > 1e10 missing since no atmospheric data should be this large.
    MissingVal = 1e10
    for i in range(len(Var)):
        if (P[i] > MissingVal    or Var[i] > MissingVal or
            VarErr[i] > MissingVal or VarQc[i] > MissingVal):
            continue

        ObsDict[cm.NC_LAT_NAME].append(ObsLat[i])
        ObsDict[cm.NC_LON_NAME].append(ObsLon[i])
        ObsDict[cm.NC_DATE_STAMP_NAME].append(AbsDate[i])
        ObsDict[cm.NC_TIME_STAMP_NAME].append(AbsTime[i])
        ObsDict[cm.NC_TIME_NAME].append(EpochTime[i])
        ObsDict[cm.NC_P_NAME].append(P[i])
        ObsDict[Vname].append(Var[i])
        ObsDict[VerrName].append(VarErr[i])
        ObsDict[VqcName].append(VarQc[i])

    return ObsDict


############################################################################
# CLASSES
############################################################################

# The BUFR format is extremely flexible, and different obs types have taken
# advantage of that fact. This has resulted in the requirement of utilizing
# different algorithms to extract obs data for different obs types. Ie, it's
# extremely difficult to force the different formats into a common algorithm.
# Using a base class with a simple extraction algorithm which can be overridden
# in a derived class seems to be a good way to handle this situation.
#
# For the extraction, it does appear that many obs types will place a header
# at the front of an BUFR subset which consists of a simple list of BUFR
# mnemonics. The header is followed by the obs data which can be a simple
# list of mnemonics, but typically is a more complex structure with 
# replications, sequences and events. The header extraction algorithm can
# (for now) belong in the base class, and the obs data extraction algorithms
# can belong in the derived classes (ie, specific to each obs type).
#
# Define the base class with a simple method that assumes all variables have a
# one-to-one corrspondence with a BUFR mnemonic. More complex examples can
# override the convert() method or its sub-methods.
#
# The format for an entry in the *_spec lists is:
#
#    [ nc_varname, mnemonic, data_type, dim_names, dim_sizes, created ] 
#
#        nc_varname: netcdf variable name
#        mnemonic:   BUFR mnemonic
#        data_type:  float, integer, string, ...
#        dim_names:  (list of dimension names)
#        dim_sizes:  (list of dimension sizes)
#        created:    flag: True  - nc variable has been created
#                          False - nc variable has not been created
#

################################# Base Observation Type ############################
class ObsType(object):
    ### Constructor ###
    def __init__(self):
        self.bufr_ftype = cm.BFILE_UNDEF
        self.mtype_re = 'UnDef'
        self.obs_type = 'UnDef'

    ### methods ###

    ######################################################################
    # This method is a dummy method that returns an error. It is intended
    # to be run when an obs type does not have its own method defined yet.
    def BufrToDict(self, *Args):
        print("ERROR: BufrToDict method is not defined for obs type: {0:s}".format(self.obs_type))
        sys.exit(1)

    ######################################################################
    # This method is a dummy method that returns an error. It is intended
    # to be run when an obs type does not have its own method defined yet.
    def PrepBufrToDict(self, *Args):
        print("ERROR: PrepBufrToDict method is not defined for obs type: {0:s}".format(self.obs_type))
        sys.exit(1)

################################# Aircraft Observation Type ############################
class AircraftObsType(ObsType):
    ### initialize data elements ###
    def __init__(self, bf_type):
        super(AircraftObsType, self).__init__()

        self.bufr_ftype = bf_type
        self.obs_type = 'Aircraft'

        if (bf_type == cm.BFILE_BUFR):
            self.mtype_re = '^NC00400[14]'
        elif (bf_type == cm.BFILE_PREPBUFR):
            self.mtype_re = 'AIRC[AF][RT]'

    ### methods ###
    
################################# Radiosonde Observation Type ############################
class SondesObsType(ObsType):
    ### initialize data elements ###
    def __init__(self, bf_type):
        super(SondesObsType, self).__init__()

        self.bufr_ftype = bf_type
        self.obs_type = 'Sondes'

        if (bf_type == cm.BFILE_BUFR):
            self.mtype_re = 'UnDef'
        elif (bf_type == cm.BFILE_PREPBUFR):
            self.mtype_re = 'ADPUPA'

    ### methods ###
    
    ###############################################################################
    # This method will read the current subset from the prepBUFR file and convert
    # the data to a table where each row is a single observation. Metadata will
    # be added in columns in the table in order to identify grouping such as 
    # channels in satellite obs, or levels in radiosonde obs, etc.
    #
    # The conversion is accomplished in two steps:
    #   1) Load BUFR data into a dictionary structure
    #   2) List out the data in the dictionary structure one obs at a time
    #
    def PrepBufrToDict(self, Bfid, MsgType, MsgDate):
        # Get the header data
        #    index   mnemonic
        #      0       SID        Station ID
        #      1       XOB        Longitude
        #      2       YOB        Latitude
        #      3       DHR        Time offset (from MsgDate)
        #
        # The Lon, Lat, Time values identify a unique sounding, air flight, etc.
        # Call this object a "record". Make a dictionary key (record) from these
        # that can be used for combining subsets belonging to the same record.
        #
        # The read_subset() method will return a masked array containing the requested
        # data. The array will be shaped: (number of mnemonics, 1).
        BufrVals = Bfid.read_subset("SID XOB YOB DHR").data
        StationId     = BufrFloatToString(BufrVals[0,0])
        HdrLon        = BufrVals[1,0]
        HdrLat        = BufrVals[2,0]
        HdrTimeOffset = BufrVals[3,0]

        # Form the record key
        RecKey = "{0:f}:{1:f}:{2:f}".format(HdrLon, HdrLat, HdrTimeOffset)

        # Get the observation data
        #    index   mnemonic
        #      0       XDR        Longitude
        #      1       YDR        Latitude
        #      2       HRDR       Time offset (from MsgDate)
        #      3       POB        Pressure
        #
        #      4       TOB        Temperature
        #      5       TOE        Temperature error
        #      6       TQM        Temperature quality code
        #
        #      7       QOB        Specific humidity
        #      8       QOE        Specific humidity error
        #      9       QQM        Specific humidity quality code
        #
        #     10       UOB        Eastward wind
        #     11       VOB        Northward wind
        #     12       WOE        Wind error (u,v combined)
        #     13       WQM        Wind quality code (u,v combined)
        #
        # The Locations data are the Lat, Lon, Time
        # values used in the DA process. These are unique for each observation (point)
        # in conventional data. A radiosonde drifts as it rises giving each measurement
        # in the sounding a unique lat, lon, time. Same thing for measurements taken
        # from an aircraft.
        BufrVals = Bfid.read_subset("XDR YDR HRDR POB TOB TOE TQM QOB QOE QQM UOB VOB WOE WQM").data
        ObsLon        = BufrVals[0,:]
        ObsLat        = BufrVals[1,:]
        ObsTimeOffset = BufrVals[2,:]
        P             = BufrVals[3,:]
        T             = BufrVals[4,:]
        Terr          = BufrVals[5,:]
        Tqc           = BufrVals[6,:]
        Q             = BufrVals[7,:]
        Qerr          = BufrVals[8,:]
        Qqc           = BufrVals[9,:]
        U             = BufrVals[10,:]
        V             = BufrVals[11,:]
        Werr          = BufrVals[12,:]
        Wqc           = BufrVals[13,:]

        # Convert the MsgDate, ObsTimeOffset to an absolute time.
        [ AbsDate, AbsTime, EpochTime ] = DateOffsetToAbsTime(MsgDate, ObsTimeOffset)

        # Form the dictionary. Organize it first by the RecKey and then by variable.
        # Treat the StationId as a variable with one value. The other variables
        # (t,u,q,v) each get a list of point data. Include Lat, Lon, Time,
        # Pressure (level), ObsValue, ObsError and ObsQcode.
        ObsDict = { RecKey : { 'VARIABLES' : {}, 'METADATA' : {} } }

        # Add in metadata
        ObsDict[RecKey]['METADATA'][cm.NC_STATION_ID_NAME] = [ StationId ]

        # Add in variables
        # Temperature
        ObsDict[RecKey]['VARIABLES'][cm.NC_T_NAME] = PrepBufrDictCreateVar(cm.NC_T_NAME,
            ObsLat, ObsLon, AbsDate, AbsTime, EpochTime, P, T, Terr, Tqc)

        # Specific humidity
        ObsDict[RecKey]['VARIABLES'][cm.NC_Q_NAME] = PrepBufrDictCreateVar(cm.NC_Q_NAME,
            ObsLat, ObsLon, AbsDate, AbsTime, EpochTime, P, Q, Qerr, Qqc)

        # Zonal wind
        ObsDict[RecKey]['VARIABLES'][cm.NC_U_NAME] = PrepBufrDictCreateVar(cm.NC_U_NAME,
            ObsLat, ObsLon, AbsDate, AbsTime, EpochTime, P, U, Werr, Wqc)
         
        # Meridional wind
        ObsDict[RecKey]['VARIABLES'][cm.NC_V_NAME] = PrepBufrDictCreateVar(cm.NC_V_NAME,
            ObsLat, ObsLon, AbsDate, AbsTime, EpochTime, P, V, Werr, Wqc)
         
        return ObsDict

########################### Radiance (AMSU-A) Observation Type ############################
class AmsuaObsType(ObsType):
    ### initialize data elements ###
    def __init__(self, bf_type):
        super(AmsuaObsType, self).__init__()

        self.bufr_ftype = bf_type
        self.obs_type = 'Amsua'

        if (bf_type == cm.BFILE_BUFR):
            self.mtype_re = '^NC021023'
        elif (bf_type == cm.BFILE_PREPBUFR):
            self.mtype_re = 'UnDef'

    ### methods ###
    
########################### GPSRO Observation Type ############################
class GpsroObsType(ObsType):
    ### initialize data elements ###
    def __init__(self, bf_type):
        super(GpsroObsType, self).__init__()

        self.bufr_ftype = bf_type
        self.obs_type = 'Gpsro'

        if (bf_type == cm.BFILE_BUFR):
            self.mtype_re = '^NC003010'
        elif (bf_type == cm.BFILE_PREPBUFR):
            self.mtype_re = 'UnDef'

    ### methods ###

