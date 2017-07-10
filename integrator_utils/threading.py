import threading


########################################
def get_thread_name():
    return threading.currentThread().getName()


# don't know how to do this if there are no other_args,
# except by passing it an empty list in th place

########################################
def parallelize(no_threads, embarassingly_pllbl_fn, list, other_args):
    if (no_threads < 1):
        print "number of threads is expected to be >= 1"
        return False

    if (no_threads == 1):
        ret = embarassingly_pllbl_fn(list, other_args)
        return ret

    #########################################
    # nontrivial
    load = []
    for thr in range(no_threads):
        load.append(0)

    for job in range(len(list)):
        load[(job % no_threads)] += 1

    # run
    total = 0
    for thr in range(no_threads):
        thr_from = total
        thr_to = total + load[thr]
        total += load[thr]

        if (thr_from >= len(list)):
            break
        if (thr == no_threads - 1):
            thr_to = len(list)

        thread = threading.Thread(target=embarassingly_pllbl_fn,
                                  args=(list[thr_from:thr_to], other_args))
        try:
            thread.start()
        except:
            print "Error: unable to start thread"
            return False



