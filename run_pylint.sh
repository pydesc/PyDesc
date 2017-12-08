#!/bin/sh
python -m pylint.lint -d I,C0301,R0901,R0902,R0903,R0904,R0913,R0915,W0141,W0142,W0232,W0613 pydesc --msg-template='{path}:{line}: [{msg_id}, {obj}] {msg}'
echo

