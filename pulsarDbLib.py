def generate(env, **kw):
    env.Tool('addLibrary', library = ['pulsarDb'], package = 'pulsarDb')
    env.Tool('st_appLib')
    env.Tool('st_facilitiesLib')
    env.Tool('st_streamLib')
    env.Tool('timeSystemLib')
    env.Tool('tipLib')

def exists(env):
    return 1
