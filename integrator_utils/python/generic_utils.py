import inspect

def frame_string():
  callerframerecord = inspect.stack()[1]    # 0 represents this line
                                            # 1 represents line at caller
  frame = callerframerecord[0]
  info = inspect.getframeinfo(frame)
  return "{}:{}():{}".format(info.filename, info.function, info.lineno)
