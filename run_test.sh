#!/usr/bin/env bash

while [[ $# -gt 1 ]]
do
_key="$1"

case $_key in
    -t|--tests)
        _tests="$2"
        shift
        ;;
    *)
        echo "Unknown option: $_key"
        exit 1
        ;;
esac
shift # past argument or value
done

if [ "x$_tests" = "x" ]; then
    _tests=`ls tests/test_*.cpp | sed -e 's/tests\/\(.*\)\.cpp/\1/'`
fi

echo $_tests

case "$_tests" in
    *unit*)
        _dc_run_opts="--no-deps --rm"
        ;;
    *)
        _dc_run_opts="--rm"
        ;;
esac

_dc_opts="-f docker-compose.yml"
exit_code=0
for _t in $_tests; do
    _command="docker-compose $_dc_opts run $_dc_run_opts hssp ./$_t"
    echo $_command
    $_command
    e=$?
    if [ $e -ne 0 ] ; then
        exit_code=$e
    fi
done

# Remove all containers and network only if the tests passed. Keep them around
# for debugging the failed tests.
if [ $exit_code -eq 0 ]; then
    docker-compose down
fi
exit $exit_code

