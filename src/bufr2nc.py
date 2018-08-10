#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import re
import argparse

import ncepbufr

import netCDF4 as nc

import bufr2ncCommon as cm
import bufr2ncObsTypes as ot

###################################################################################
# CLASSES
###################################################################################

################################################################################
# MessageSelector
#
# This class is used to apply the BUFR message selection criteria while reading
# a BUFR file.
class MessageSelector(object):

    ###############################################################################
    # Contstructor
    def __init__(self, MtypeRe, MaxNumMsg, ThinInterval):
        self.mtype_re = MtypeRe
        self.max_num_msg = MaxNumMsg
        self.thin_interval = ThinInterval

        self.num_msg_selected = 0
        self.num_msg_mtype = 0

    ###############################################################################
    # This method is the message selector. It will apply selection filters
    # to the input BUFR messages.
    def select_next_msg(self, Bfid):
        got_a_msg = False
        # Grab the next message
        while (Bfid.advance() == 0):
            # Skip this message if not the desired type
            if (re.search(self.mtype_re, Bfid.msg_type)):
                # Keep count of the messages that match the desired type, which is
                # needed to do the selection filtering.
                self.num_msg_mtype += 1

                # Apply the filtering. Default is to take all messages
                Select = True

                # If the max_num_msg parameter is greater than zero, then use it to limit
                # the number of messages that are selected.
                if (self.max_num_msg > 0):
                    Select = (self.num_msg_selected < self.max_num_msg)

                # If the thinning interval is greater than 1, then use it to further select
                # every n-th message.
                if (self.thin_interval > 1):
                    Select = Select and ((self.num_msg_mtype % self.thin_interval) == 0)

                # If Select is true, the current message has been selected. Keep
                # track of how many messages have been selected, plus break out of
                # the loop and return.
                if (Select):
                    self.num_msg_selected += 1
                    got_a_msg = True
                    break

        return got_a_msg
 

###################################################################################
# SUBROUTINES
###################################################################################

#############################################################################
# This routine will read through the BUFR file, select messages according to
# the input selection parameters, convert the BUFR messages from a list
# of subsets to a table and write the tables into the netcdf file.
def ConvertMessagesToTables(BufrFname, NetcdfFname, Obs, MaxNumMsg, ThinInterval):
    # Create a selector object
    Msel = MessageSelector(Obs.mtype_re, MaxNumMsg, ThinInterval)

    # Open files
    Nfid = nc.Dataset(NetcdfFname, 'w', format='NETCDF4')
    Bfid = ncepbufr.open(BufrFname)

    # Do the message selection according to the parameters:
    #   MaxNumMsg - maximum number of messages (if <0, then ignore this spec)
    #   ThinInterval - select every n-th message

    NumMsgs = 0
    NumSubsets = 0

    ObsDict = { }
    RecSet = set()
    MdataSet = set()
    VarSet = set()

    while(Msel.select_next_msg(Bfid)):
        NumMsgs += 1

        MsgType = Bfid.msg_type
        MsgDate = Bfid.msg_date

        # Walk though all subsets, collect obs data and store in a dictionary.
        # The dictionary is used so that multiple subsets with the same Lat, Lon, Time
        # signature can be combined into one record.
        while (Bfid.load_subset() == 0):
            NumSubsets += 1

            # Create a dictionary for this subset holding the observation values
            # extracted from the BUFR subset.
            if (Obs.bufr_ftype == cm.BFILE_PREPBUFR):
                SubsetObsDict = Obs.PrepBufrToDict(Bfid, MsgType, MsgDate)
            elif (Obs.bufr_ftype == cm.BFILE_BUFR):
                SubsetObsDict = Obs.BufrToDict(Bfid, MsgType, MsgDate)

            # Append the subset dictionary to the message dictionary
            for RecKey in SubsetObsDict:
                RecSet.add(RecKey)
                if RecKey not in ObsDict:
                    ObsDict[RecKey] = { 'METADATA' : { }, 'VARIABLES' : { } }

                # Add in metadata values
                for MdataKey in SubsetObsDict[RecKey]['METADATA']:
                    MdataSet.add(MdataKey)
                    if MdataKey not in ObsDict[RecKey]['METADATA']:
                        # Add in subset metadata value (list)
                        ObsDict[RecKey]['METADATA'][MdataKey] = SubsetObsDict[RecKey]['METADATA'][MdataKey]
                    else:
                        # Append subset metadata value (list)
                        ObsDict[RecKey]['METADATA'][MdataKey] += SubsetObsDict[RecKey]['METADATA'][MdataKey]

                # Add in variable values
                for VarKey in SubsetObsDict[RecKey]['VARIABLES']:
                    VarSet.add(VarKey)
                    if VarKey not in ObsDict[RecKey]['VARIABLES']:
                        # Add in subset variable values (dict containing lists)
                        ObsDict[RecKey]['VARIABLES'][VarKey] = SubsetObsDict[RecKey]['VARIABLES'][VarKey]
                    else:
                        # Append subset variable values (dict containing lists)
                        for ObsKey in SubsetObsDict[RecKey]['VARIABLES'][VarKey]:
                            ObsDict[RecKey]['VARIABLES'][VarKey][ObsKey] += SubsetObsDict[RecKey]['VARIABLES'][VarKey][ObsKey]


    # Convert the sets of keys to sorted lists. The numbering for the table
    # will be taken from these lists.
    RecNames   = sorted(list(RecSet))
    MdataNames = sorted(list(MdataSet))
    VarNames   = sorted(list(VarSet))
    
    # Write the data out into the netcdf file. Use groups to organize the data
    # in the file according to records, metadata and variables.

    for irec in range(len(RecNames)):
        RecNum = irec + 1
        RecKey = RecNames[irec]
        RecName = "Rec_" + "{0:d}".format(RecNum)

        Rgroup = Nfid.createGroup(RecName)

        MdataDict = ObsDict[RecKey]['METADATA']
        VarDict   = ObsDict[RecKey]['VARIABLES']

        for imdata in range(len(MdataNames)):
            MdataKey = MdataNames[imdata]
            MdataVal = MdataDict[MdataKey]

            Rgroup.setncattr(MdataKey, MdataVal)

        for ivar in range(len(VarNames)):
            VarNum = ivar + 1
            VarKey = VarNames[ivar]

            VarGroup = Rgroup.createGroup(VarKey)

            Nvals = len(VarDict[VarKey][VarKey])
            VarGroup.createDimension("nvals", Nvals)

            NcVar = VarGroup.createVariable(cm.NC_LAT_NAME, "f4", ("nvals"))
            NcVar[:] = VarDict[VarKey][cm.NC_LAT_NAME]

            NcVar = VarGroup.createVariable(cm.NC_LON_NAME, "f4", ("nvals"))
            NcVar[:] = VarDict[VarKey][cm.NC_LON_NAME]

            NcVar = VarGroup.createVariable(cm.NC_DATE_STAMP_NAME, "i4", ("nvals"))
            NcVar[:] = VarDict[VarKey][cm.NC_DATE_STAMP_NAME]

            NcVar = VarGroup.createVariable(cm.NC_TIME_STAMP_NAME, "i4", ("nvals"))
            NcVar[:] = VarDict[VarKey][cm.NC_TIME_STAMP_NAME]

            NcVar = VarGroup.createVariable(cm.NC_TIME_NAME, "f8", ("nvals"))
            NcVar[:] = VarDict[VarKey][cm.NC_TIME_NAME]

            NcVar = VarGroup.createVariable(cm.NC_P_NAME, "f4", ("nvals"))
            NcVar[:] = VarDict[VarKey][cm.NC_P_NAME]

            NcVar = VarGroup.createVariable(VarKey, "f4", ("nvals"))
            NcVar[:] = VarDict[VarKey][VarKey]

            NcVar = VarGroup.createVariable(VarKey+'_err', "f4", ("nvals"))
            NcVar[:] = VarDict[VarKey][VarKey+'_err']

            NcVar = VarGroup.createVariable(VarKey+'_qc', "f4", ("nvals"))
            NcVar[:] = VarDict[VarKey][VarKey+'_qc']


    # Finish up
    Bfid.close()

    Nfid.sync()
    Nfid.close()

    # return counts
    return [ NumMsgs, NumSubsets ]

###################################################################################
# MAIN
###################################################################################
ScriptName = os.path.basename(sys.argv[0])

# Parse command line
ap = argparse.ArgumentParser()
ap.add_argument("obs_type", help="observation type")
ap.add_argument("input_bufr", help="path to input BUFR file")
ap.add_argument("output_netcdf", help="path to output netCDF4 file")
ap.add_argument("-m", "--maxmsgs", type=int, default=-1,
                help="maximum number of messages to keep", metavar="<max_num_msgs>")
ap.add_argument("-t", "--thin", type=int, default=1,
                help="select every nth message (thinning)", metavar="<thin_interval>")
ap.add_argument("-c", "--clobber", action="store_true",
                help="allow overwrite of output netcdf file")
ap.add_argument("-p", "--prepbufr", action="store_true",
                help="input BUFR file is in prepBUFR format")

MyArgs = ap.parse_args()

ObsType = MyArgs.obs_type
BufrFname = MyArgs.input_bufr
NetcdfFname = MyArgs.output_netcdf
MaxNumMsg = MyArgs.maxmsgs
ThinInterval = MyArgs.thin
ClobberOfile = MyArgs.clobber
if (MyArgs.prepbufr):
    BfileType = cm.BFILE_PREPBUFR
else:
    BfileType = cm.BFILE_BUFR

# Check files
BadArgs = False
if (not os.path.isfile(BufrFname)): 
    print("ERROR: {0:s}: Specified input BUFR file does not exist: {1:s}".format(ScriptName, BufrFname))
    print("")
    BadArgs = True

if (os.path.isfile(NetcdfFname)):
    if (ClobberOfile):
        print("WARNING: {0:s}: Overwriting nc file: {1:s}".format(ScriptName, NetcdfFname))
        print("")
    else:
        print("ERROR: {0:s}: Specified nc file already exists: {1:s}".format(ScriptName, NetcdfFname))
        print("ERROR: {0:s}:   Use -c option to overwrite.".format(ScriptName))
        print("")
        BadArgs = True

# Check for valid observation type, and if okay instantiate an obs type object.
if (ObsType == 'Aircraft'):
    Obs = ot.AircraftObsType(BfileType)
elif (ObsType == 'Sondes'):
    Obs = ot.SondesObsType(BfileType)
elif (ObsType == 'Amsua'):
    Obs = ot.AmsuaObsType(BfileType)
elif (ObsType == 'Gpsro'):
    Obs = ot.GpsroObsType(BfileType)
else:
    print("ERROR: {0:s}: Unknown observation type: {1:s}".format(ScriptName, ObsType))
    print("")
    BadArgs = True

if (not BadArgs):
    if (Obs.mtype_re == 'UnDef'):
        if (BfileType == cm.BFILE_BUFR):
            print("ERROR: {0:s}: Observation type {1:s} for BUFR format is undefined".format(ScriptName, ObsType))
        elif (BfileType == cm.BFILE_PREPBUFR):
            print("ERROR: {0:s}: Observation type {1:s} for prepBUFR format is undefined".format(ScriptName, ObsType))
        print("")
        BadArgs = True

if (BadArgs):
    sys.exit(2)

# Arguments are okay
print("Converting BUFR to netCDF")
print("  Observation Type: {0:s}".format(ObsType))
if (BfileType == cm.BFILE_BUFR):
    print("  Input BUFR file (BUFR format): {0:s}".format(BufrFname))
elif (BfileType == cm.BFILE_PREPBUFR):
    print("  Input BUFR file (prepBUFR format): {0:s}".format(BufrFname))
print("  Output netCDF file: {0:s}".format(NetcdfFname))
if (MaxNumMsg > 0):
    print("  Limiting nubmer of messages to record to {0:d} messages".format(MaxNumMsg))
if (ThinInterval > 1):
    print("  Thining: selecting every {0:d}-th message".format(ThinInterval))
print("")

# All BUFR files are organized as a list of messages. Under each message is a list of
# subsets. Each subset holds data for an observation set, such as a sounding for sondes or
# the brightness temperatures for all channels in a satellite sensor. The message, subset
# organization is the same in all BUFR files, but the format for a subset changes for
# each observation type. So the idea here is to walk through messages and subsets in the
# same way regardless of observation type, and use a class structure (by obs type) for the
# parsing of the subset.
#
# BUFR messages can contain related data, such as two subsets that combine to form one
# radiosonde sounding. For this reason, use messages as the smallest unit to be worked
# upon. For example, dole out messages to multiple processes in an MPI run.
#
# Break the BUFR to netcdf conversion into two steps:
#  1) convert each BUFR message to a table and store the tables in the netcdf file
#  2) run through the netcdf file and combine the tables into one larger table

print("Converting BUFR messages to obs tables: ")
[ NumMsg, NumSubsets ] = ConvertMessagesToTables(BufrFname, NetcdfFname,
                                                 Obs, MaxNumMsg, ThinInterval)
print("  Number of messages selected: {0:d}".format(NumMsg))
print("  Number of subsets converted: {0:d}".format(NumSubsets))
print("")



sys.exit(0)
