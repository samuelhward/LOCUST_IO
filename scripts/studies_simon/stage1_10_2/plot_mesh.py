#change path to the relevant study output folder
# use #949494 as below colour

import pathlib

simulations=[
    '*-3_p*90000*',
    #'*-4_p*90000*',
    #'*-3_-6*90000*',
    #'*-4_-5*90000*',
    #'*-3_p*60000*',
    #'*-4_p*60000*',
    #'*-3_-6*60000*',
    #'*-4_-5*60000*',
    #'*phaseu_20.033333333333328_phasem_86.7_phasel_40.033333333333*',
    #'*phaseu_-19.966666666666672_phasem_46.7_phasel_0.03333333333333588*',
    ] # directory search patterns
components=[
            '01', #first wall
            '02', #dome
            '03', #inner vertical supports
            '04', #outer vertical supports
            '05', #horizontal supports
            '06', #inner plate
            '07', #outer plate
            '08', #outer pipes
            '09', #inner pipes
            '10', #outer strikepoint
            '11', #inner strikepoint
            '12', #base
                ] # component IDs

""" view numbers
1 - look at PFC wall loads
"""
view_number=-1
camera_views={}
camera_views['position']=[
    (-80.6175711892685,-25.17004340407864,7.653252517525814),
]
camera_views['view_up']=[
    (0,0,1)
]
camera_views['view_angle']=[
    10
]
camera_views['focal_point']=[
    (40.323941219165334, 14.862391198241832, 0.06896562340792889),
]
camera_views['azimuth']=[
    0,
]
camera_views['elevation']=[
    0,
]
camera_views['roll']=[
    0,
]

def set_camera(views,position,view_up,view_angle,focal_point,azimuth,elevation,roll):
    for view in views:
        SetActiveView(view)
        camera=view.GetActiveCamera() 
        camera.SetPosition(*position) #GetPosition() 
        camera.SetViewUp(*view_up)  
        camera.SetViewAngle(view_angle) #GetViewAngle()
        camera.SetFocalPoint(focal_point) #GetFocalPoint()
        camera.Azimuth(azimuth)
        camera.Elevation(elevation)
        camera.Roll(roll)
        Render()

def render_mesh(simulations,components,vminmax=[0,0.5],opacity_transfer=False,add_label=True,colour_by='Power_Flux_[MW/m2]'):
    import datetime
    if not components:
        patterns=['*_1.vtk']
    else:
        patterns=[f'*_00{component}.vtk' for component in components]
    # paths to output directories assumed cwd
    # set to remove degeneracy
    dir_names={path for simulation in simulations for path in pathlib.Path('.').glob(simulation)} 
    dirs={}
    for dir_name in dir_names:
        files=[]
        for pattern in patterns:
            #find earliest file as likely contains the real PFC power loss data
            date=datetime.datetime.today()
            component_files=pathlib.Path(dir_name).glob(pattern)
            target_filename=''
            for component_file in component_files:
                date_file=datetime.datetime.strptime(str(component_file.parts[-1]).split('_')[2],'%d-%m-%Y')
                #print(date_file)
                if date_file<date: 
                    date=date_file
                    target_filename=component_file
            files.append(target_filename)
        dirs[dir_name] = set(files)
    layout=GetLayout()
    views=GetRenderViews()
    for _ in range(len(dirs)-1):
        view=CreateRenderView()
        views.append(view)
        AssignViewToLayout(view=view, layout=layout, hint=_+1)
    labels=[]
    for counter,((dir,filepaths),view) in enumerate(zip(dirs.items(),views)): 
        print(f'plotting {dir}...')
        SetActiveView(view)
        if add_label:
            label=Text()
            label_text=r'$n$='+'+'.join(str(abs(int(mode))) for mode in str(dir).split('ntor')[1].split('phaseu')[0].split('_') if mode)
            label_text+=r', $\Phi_{\mathrm{u,m,l}}=$'
            label_text+=','.join(str(int(float(str(dir).split(f'phase{row}')[-1].split('_')[1]))) for row in ['u','m','l'])
            label.Text=label_text
            labels.append(label)
            Show(label,view=view)
        for filepath in filepaths:
            print(f'plotting {filepath.parts[-1]}...')
            reader=OpenDataFile(str(filepath))
            UpdatePipeline()
            display=Show(reader,view=view)
            ColorBy(display,('CELLS',colour_by))
        Render(view=view)
        colorMap = GetColorTransferFunction(colour_by)
        if vminmax:
            colorMap.RescaleTransferFunction(*[minmax for minmax in vminmax])
        if opacity_transfer: 
            colorMap.EnableOpacityMapping = 1
            opacityMap = GetOpacityTransferFunction(colour_by)
            opacityMap.RescaleTransferFunction(0.001, 1)
            opacityMap.Points=[
                0.,0.000,0.5,0.,
                0.01,1.,0.5,0.,
                1.,1.,0.5,0.,
            ]

render_mesh(simulations,components,[0.001,0.1],opacity_transfer=True)

render_mesh(simulations,components,opacity_transfer=False,add_label=False,colour_by='Component_ID',vminmax=[0,12])

#set_camera(views=GetRenderViews(),**{key:value[0] for key,value in camera_views.items()})
#colorMap = GetColorTransferFunction('Power_Flux_[MW/m2]')
colorMap = GetColorTransferFunction('Power_Flux_[MW/m2]')
colorMap.RescaleTransferFunction(0.001,0.3)
colorMap.EnableOpacityMapping = 0
opacityMap = GetOpacityTransferFunction('Power_Flux_[MW/m2]')
opacityMap.RescaleTransferFunction(0.001, 1)
opacityMap.Points=[
    0.,0.000,0.5,0.,
    0.001,1.,0.5,0.,
    1.,1.,0.5,0.,
]

#display = GetDisplayProperties(display, view)
#display.SetScalarBarVisibility(view, True)

layout=GetLayout()
SaveScreenshot("paper_2_CAD_model_key.png",ImageResolution=[5000, 2500],layout=layout)




