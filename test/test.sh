#! /usr/bin/env bash

function header()   { echo -e "\n\033[1m$@\033[0m"; }
function success()  { echo -e " \033[1;32m==>\033[0m  $@"; }
function error()    { echo -e " \033[1;31mX\033[0m  $@"; }
function arrow()    { echo -e " \033[1;34m==>\033[0m  $@"; }

TEST_NAME=doTest.sh
GLOBALS_FILE=globals.conf

#DEFAULT VALUES
RUN_COMMAND=mpirun
TEST_CLASS=essential
CONFIG=icc
CC4S_PATH=


__SCRIPT_VERSION="0.1"
__SCRIPT_NAME=$( basename $0 )
__DESCRIPTION="Test suite for cc4s"
__OPTIONS=":hvt:r:c:x:"

function usage_head() { echo "Usage :  $__SCRIPT_NAME [-h|-help] [-v|-version] [-c CLASS_NAME] [-r RUN_COMMAND]"; }
function usage ()
{
cat <<EOF
$(usage_head)

    $__DESCRIPTION

    Options:
    -h|help       Display this message
    -v|version    Display script version
    -r            Run command
                  (e.g. "mpirun", "mpiexec.hydra -bootstrap ll" ...)
                  Default: ${RUN_COMMAND}
    -t            Test class (e.g. all essential extended)
                  Default: ${TEST_CLASS}
    -c            Compiled configuration for the cc4s binary
                  Default: ${CONFIG}
    -x            Path to cc4s executable, this overrides -c flag

EOF
}    # ----------  end of function usage  ----------

while getopts $__OPTIONS opt
do
  case $opt in

  h|help     )  usage; exit 0   ;;

  v|version  )  echo "$__SCRIPT_NAME -- Version $__SCRIPT_VERSION"; exit 0   ;;

  r  )
    RUN_COMMAND=${OPTARG} 
    arrow "Using RUN_COMMAND = ${RUN_COMMAND}"
    ;;

  c  )
    CONFIG=${OPTARG}
    arrow "Using CONFIG = ${CONFIG}"
    ;;

  c  )
    CC4S_PATH=${OPTARG}
    arrow "Using CC4S_PATH = ${CC4S_PATH}"
    ;;

  t  )
    TEST_CLASS=${OPTARG}
    arrow "Using TEST_CLASS = ${TEST_CLASS}"
    ;;

  * )  echo -e "\n  Option does not exist : $OPTARG\n"
      usage_head; exit 1   ;;

  esac    # --- end of case ---
done
shift $(($OPTIND-1))

# Check if cc4s path was overriden
if [[ -z ${CC4S_PATH} ]] ; then
  CC4S_PATH=$(readlink -f "../build/${CONFIG}/bin/Cc4s")
  arrow "Using CC4S_PATH = ${CC4S_PATH}"
fi

# Check for cc4s executable
if [[ ! -x ${CC4S_PATH} ]] ; then
  error "CC4S not found"
  exit 1
fi

#if the class is all, then all test cases should
#match
if [[ ${TEST_CLASS} = "all" ]]; then
  TEST_CLASS="@CLASS"
fi

arrow "Sourcing ${GLOBALS_FILE} file"
source ${GLOBALS_FILE}

MAIN_TEST_FOLDER=$(pwd)
TEST_OUT_FILE=test.out
FAILED_TEST_COUNT=0
TEST_COUNT=0

for TEST_SCRIPT in $(find . -name ${TEST_NAME}); do
  if grep "@CLASS" ${TEST_SCRIPT} | grep ${TEST_CLASS} &> /dev/null; then
    let TEST_COUNT+=1
  else
    continue
  fi
  TEST_RESULT=1
  header "Testing ${TEST_SCRIPT} ... "
  TEST_FOLDER=$(dirname ${TEST_SCRIPT})
  cd ${TEST_FOLDER}
  source ${TEST_NAME} > ${TEST_OUT_FILE}
  if [[ ${TEST_RESULT} = 0 ]]; then
    success "Sucess"
  else
    error "Test FAILED"
    let FAILED_TEST_COUNT=+1
  fi
  cd ${MAIN_TEST_FOLDER}
done

header "${TEST_COUNT} tests DONE for class '${TEST_CLASS}'"
header "${FAILED_TEST_COUNT} tests FAILED for class '${TEST_CLASS}'"

if [[ ${FAILED_TEST_COUNT} != 0 && ${TEST_CLASS} == "essential" ]]; then
  cat <<EOF

############################################
#  DO NOT COMMIT TO MASTER BRANCH, PLEASE  #
############################################

EOF
fi


