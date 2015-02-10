
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 4.3.1 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 4.3.1

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1133, 550]
      renderView1.InteractionMode = '2D'
      renderView1.OrientationAxesVisibility = 0
      renderView1.CenterOfRotation = [0.630069932613178, 3.4715622348802397, 0.0]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [0.630069932613178, 3.4715622348802397, 10000.0]
      renderView1.CameraFocalPoint = [0.630069932613178, 3.4715622348802397, 0.0]
      renderView1.CameraParallelScale = 88.80059596465395
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='OW1_global_%t.png', freq=1, fittoscreen=0, magnification=1, width=1133, height=550)

      # Create a new 'Render View'
      renderView2 = CreateView('RenderView')
      renderView2.ViewSize = [1133, 550]
      renderView2.OrientationAxesVisibility = 0
      renderView2.CenterOfRotation = [32.42239624147474, -32.36633541180826, 0.0]
      renderView2.StereoType = 0
      renderView2.CameraPosition = [32.42239624147474, -33.28006052384497, 62.51792327499739]
      renderView2.CameraFocalPoint = [32.42239624147474, -33.28006052384497, 0.0]
      renderView2.CameraParallelScale = 265.6939056222687
      renderView2.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView2,
          filename='OW2_Agulhas_%t.png', freq=1, fittoscreen=0, magnification=1, width=1133, height=550)

      # Create a new 'Render View'
      renderView3 = CreateView('RenderView')
      renderView3.ViewSize = [1133, 550]
      renderView3.OrientationAxesVisibility = 0
      renderView3.CenterOfRotation = [-42.83581779503237, 28.038283749072495, 0.0]
      renderView3.StereoType = 0
      renderView3.CameraPosition = [-45.677635317395854, 35.71887164735211, 78.8269467380402]
      renderView3.CameraFocalPoint = [-45.677635317395854, 35.71887164735211, 0.0]
      renderView3.CameraParallelScale = 265.6939056222687
      renderView3.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView3,
          filename='OW3_Atlantic_%t.png', freq=1, fittoscreen=0, magnification=1, width=1133, height=550)

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML Partitioned Unstructured Grid Reader'
      # create a producer from a simulation input
      mpas_data_1pvtu = coprocessor.CreateProducer(datadescription, 'LON_LAT_1LAYER-primal')

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'okuboWeiss'
      okuboWeissLUT = GetColorTransferFunction('okuboWeiss')
      okuboWeissLUT.RGBPoints = [-1.0, 0.231373, 0.298039, 0.752941, 0.0, 0.865003, 0.865003, 0.865003, 1.0, 0.705882, 0.0156863, 0.14902]
      okuboWeissLUT.LockScalarRange = 1
      okuboWeissLUT.ScalarRangeInitialized = 1.0
      okuboWeissLUT.VectorMode = 'Component'

      # get opacity transfer function/opacity map for 'okuboWeiss'
      okuboWeissPWF = GetOpacityTransferFunction('okuboWeiss')
      okuboWeissPWF.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
      okuboWeissPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from mpas_data_1pvtu
      mpas_data_1pvtuDisplay = Show(mpas_data_1pvtu, renderView1)
      # trace defaults for the display properties.
      mpas_data_1pvtuDisplay.ColorArrayName = ['CELLS', 'okuboWeiss']
      mpas_data_1pvtuDisplay.LookupTable = okuboWeissLUT
      mpas_data_1pvtuDisplay.ScalarOpacityUnitDistance = 15.50910616445539

      # show color legend
      mpas_data_1pvtuDisplay.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for okuboWeissLUT in view renderView1
      okuboWeissLUTColorBar = GetScalarBar(okuboWeissLUT, renderView1)
      okuboWeissLUTColorBar.Position = [0.8429328621908128, 0.521766848816029]
      okuboWeissLUTColorBar.Position2 = [0.12, 0.43000000000000005]
      okuboWeissLUTColorBar.Title = 'okuboWeiss'
      okuboWeissLUTColorBar.ComponentTitle = '0'

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView2'
      # ----------------------------------------------------------------

      # show data from mpas_data_1pvtu
      mpas_data_1pvtuDisplay_1 = Show(mpas_data_1pvtu, renderView2)
      # trace defaults for the display properties.
      mpas_data_1pvtuDisplay_1.ColorArrayName = ['CELLS', 'okuboWeiss']
      mpas_data_1pvtuDisplay_1.LookupTable = okuboWeissLUT
      mpas_data_1pvtuDisplay_1.ScalarOpacityUnitDistance = 15.50910616445539

      # show color legend
      mpas_data_1pvtuDisplay_1.SetScalarBarVisibility(renderView2, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for okuboWeissLUT in view renderView2
      okuboWeissLUTColorBar_1 = GetScalarBar(okuboWeissLUT, renderView2)
      okuboWeissLUTColorBar_1.Title = 'okuboWeiss'
      okuboWeissLUTColorBar_1.ComponentTitle = '0'

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView3'
      # ----------------------------------------------------------------

      # show data from mpas_data_1pvtu
      mpas_data_1pvtuDisplay_2 = Show(mpas_data_1pvtu, renderView3)
      # trace defaults for the display properties.
      mpas_data_1pvtuDisplay_2.ColorArrayName = ['CELLS', 'okuboWeiss']
      mpas_data_1pvtuDisplay_2.LookupTable = okuboWeissLUT
      mpas_data_1pvtuDisplay_2.ScalarOpacityUnitDistance = 15.50910616445539

      # show color legend
      mpas_data_1pvtuDisplay_2.SetScalarBarVisibility(renderView3, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for okuboWeissLUT in view renderView3
      okuboWeissLUTColorBar_2 = GetScalarBar(okuboWeissLUT, renderView3)
      okuboWeissLUTColorBar_2.Title = 'okuboWeiss'
      okuboWeissLUTColorBar_2.ComponentTitle = '0'
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'LON_LAT_1LAYER-primal': [1, 1, 1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(False, 1)


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
