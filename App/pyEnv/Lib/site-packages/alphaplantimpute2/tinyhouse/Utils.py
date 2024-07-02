import datetime
def time_func(text):
    # This creates a decorator with "text" set to "text"
    def timer_dec(func):
        # This is the returned, modified, function
        def timer(*args, **kwargs):
            start_time = datetime.datetime.now()
            values = func(*args, **kwargs)
            print(f"{text}: {(datetime.datetime.now() - start_time).total_seconds():.2f} seconds")
            return values
        return timer

    return timer_dec
