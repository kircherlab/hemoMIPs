#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V

if [ -z "$SGE_TASK_ID" ] ; then
    echo "WARNING: NO $SGE_TASK_ID"
    exit 1
fi

CMD_FILE=$1

COMMAND=$(head -n $SGE_TASK_ID $CMD_FILE | tail -n 1)

if [ -n "$COMMAND" ] ; then
    $COMMAND
else
    echo "command \"$COMMAND\"not found"
fi

# $Id: generic_array_job.csh,v 1.1 2006/05/16 14:11:32 sbender Exp $
