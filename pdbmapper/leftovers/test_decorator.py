from decorator import tags
import sys
import time

text_start = "Start"
text_succeed = "Great!"
text_fail = "failed"
emoji = "🦸"

@tags(text_start, text_succeed, text_fail, emoji)

def second_wrap(eo):
    if int(eo) > 1: 
        time.sleep(1)
        return "hello " + eo 
    else: 
        time.sleep(1)
        raise IOError("esto no mola")

eo = sys.argv[1]
print(second_wrap(eo))