from astropy import constants as const

import numpy as np
import sys

import casac
#from taskinit import gentools
#tb = gentools(['tb'])[0]
#ms = gentools(['ms'])[0]

ms = casac.casac.ms()
tb = casac.casac.table()


####################################################
# Save visibility UV coordinates and frequencies
####################################################
def uvsaver(visname,npzname,spws,fields):
  npzname = npzname.replace('.npz','')
  if ('field-fid' not in npzname): npzname = '{0}.field-fid'.format(npzname)
  if   ('spw-sid' not in npzname): npzname = '{0}.spw-sid'.format(npzname)
  
  print(visname)
  ms.open(visname,nomodify=True)
  for spw in spws:
    for f, field in enumerate(fields):
      print('* Field {0} spw {1}'.format(field,spw))
      ms.selectinit(datadescid=spw)
      ms.select({'field_id': field})

      u = np.copy(ms.getdata(['u'])['u'])
      v = np.copy(ms.getdata(['v'])['v'])
      freqs = ms.range('chan_freq')['chan_freq'][:,0]

      uvdata = np.array([np.array([u*freqs[f]/const.c.value for f in range(len(freqs))]).flatten(),
                         np.array([v*freqs[f]/const.c.value for f in range(len(freqs))]).flatten(),
                         np.array([np.zeros(np.shape(u))+freqs[f] for f in range(len(freqs))]).flatten()])

      uvfile = npzname.replace('-fid','-{0}'.format(field)).replace('-sid','-{0}'.format(spw))+'.npz'
      np.savez_compressed(uvfile,uvdata)
      del u, v, uvdata, freqs
      ms.reset()
  ms.close()

####################################################
# Load model into a measurement set
####################################################
def uvloader(visname,npzname,todo,spws,fields):

  npzname = npzname.replace('.npz','')
  if ('field-fid' not in npzname): npzname = '{0}.field-fid'.format(npzname)
  #if   ('spw-sid' not in npzname): npzname = '{0}.spw-sid'.format(npzname)

  #tb.open('{0}/SPECTRAL_WINDOW'.format(visname),nomodify=True)
  #visfreqs = tb.getcol('CHAN_FREQ')
  #tb.close()

  visfreqs = []
  ms.open('{0}'.format(visname),nomodify=True)
  for s, spw in enumerate(spws):
    ms.selectinit(datadescid=spw)
    visfreqs.append(ms.range('chan_freq')['chan_freq'])
    ms.reset()
  ms.close()

  tb.open(visname,nomodify=False)
  visdescid = tb.getcol('DATA_DESC_ID')
  visfields = tb.getcol('FIELD_ID')

  visdata = tb.getcol('DATA')

  for s, spw in enumerate(spws):
    for f, field in enumerate(fields):
      print('* Field {0} spw {1}'.format(field,spw))
      spwindx = np.logical_and(visdescid==spw,visfields==field)

      spwname = npzname.replace('-fid','-{0}'.format(field)).replace('-sid','-{0}'.format(spw))+'.npz'
      spwload = np.load(spwname,fix_imports=True,encoding='bytes')
      spwreal = np.copy(spwload[spwload.files[0]][0].flatten()); spwreal = spwreal.reshape((visfreqs[s].shape[0],spwreal.shape[-1]/visfreqs[s].shape[0]))
      spwimag = np.copy(spwload[spwload.files[0]][1].flatten()); spwimag = spwimag.reshape((visfreqs[s].shape[0],spwimag.shape[-1]/visfreqs[s].shape[0]))
      spwload.close(); del spwload, spwname

      if   (todo== 'replace'): visdata[:,:,spwindx] = np.array([(spwreal+1.0j*spwimag),(spwreal+1.0j*spwimag)])
      elif (todo=='subtract'): visdata[:,:,spwindx] = visdata[:,:,spwindx]-np.array([(spwreal+1.0j*spwimag),(spwreal+1.0j*spwimag)])
      del spwindx, spwreal, spwimag

  tb.putcol('DATA',visdata)
  tb.flush()

  del visdescid, visfields
  del visfreqs, visdata
  tb.close()


####################################################
# Load model into a measurement set
####################################################
def uvloader2(visname,npzname,todo,spws,fields):

  npzname = npzname.replace('.npz','')
  if ('field-fid' not in npzname): npzname = '{0}.field-fid'.format(npzname)
  #if   ('spw-sid' not in npzname): npzname = '{0}.spw-sid'.format(npzname)

  ms.open(visname,nomodify=False)
  for spw in spws:
    for f, field in enumerate(fields):
      print('* Field {0} spw {1}'.format(field,spw))
      ms.selectinit(datadescid=spw)
      ms.select({'field_id': field})

      rec = ms.getdata(['data'])
      freqs = ms.range('chan_freq')['chan_freq'][:,0]

      uvfile = npzname.replace('-fid','-{0}'.format(field)).replace('-sid','-{0}'.format(spw))+'.npz'
      uvload = np.load(uvfile,fix_imports=True,encoding='bytes')

      uvreal = np.copy(uvload[uvload.files[0]][0].flatten()).reshape((len(freqs),len(rec['data'][0][0])))
      uvimag = np.copy(uvload[uvload.files[0]][1].flatten()).reshape((len(freqs),len(rec['data'][0][0])))
      uvload.close(); del uvload, uvfile

      if   (todo== 'replace'): rec['data'] = np.array([(uvreal+1.0j*uvimag),(uvreal+1.0j*uvimag)])
      elif (todo=='subtract'): rec['data'] = np.copy(rec['data'])-np.array([(uvreal+1.0j*uvimag),(uvreal+1.0j*uvimag)])
      ms.putdata(rec)
      ms.reset()
  ms.close()



