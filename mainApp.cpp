#include "mainApp.h"
#include <wx/msgdlg.h> 
#include <wx/textdlg.h> 
#include <string>

wxIMPLEMENT_APP_CONSOLE(RaytracerApp);

#include "fluid.h"

bool RaytracerApp::OnInit()
{
	int arc = wxApp::argc;
	//wxCmdLineArgsArray argv(wxApp::argv);
	for (int i=0; i<64; i++)
		engine[i].seed(i*100+1);

	/*Fluid fl(BBox(Vector(-10, -10, -10), Vector(10, 10, 10)), 64, 64, 64, 200000, 1000.);
	fl.run();
	exit(0);*/

	if (argc > 1) {
		Raytracer raytracer;
		raytracer.loadScene();
		if (argc>3)
			raytracer.load_scene(wxApp::argv[1], wxApp::argv[3]);
		else
			raytracer.load_scene(wxApp::argv[1]);
		raytracer.render_image_nopreviz();
		if (arc>2)
			save_image(wxApp::argv[2], &raytracer.image[0], raytracer.W, raytracer.H);
		exit(0);
	}


	if (!wxApp::OnInit())
		return false;

	// create the main frame window
	RaytracerFrame* frame = new RaytracerFrame();
	frame->SetSize(1200, 1000);

	renderPanel = new RenderPanel(frame, this);
	renderPanel->displayW = 1000;
	renderPanel->displayH = 800;

	wxBoxSizer * statusbar_sizer = new wxBoxSizer(wxHORIZONTAL);
	progressBar = new wxGauge(frame->GetStatusBar(), wxID_ANY, 1000, wxDefaultPosition, wxDefaultSize, wxGA_HORIZONTAL);
	infoModel = new wxStaticText(frame->GetStatusBar(), wxID_ANY, "triangles: 0, vertices: 0, bvh leaves size: 0, mouse distance: 0.000");
	infoPerf = new wxStaticText(frame->GetStatusBar(), wxID_ANY, "Time per ray: 00.000 s                                           ");
	statusbar_sizer->Add(progressBar);
	statusbar_sizer->Add(infoPerf);
	statusbar_sizer->Add(infoModel);
	frame->GetStatusBar()->SetSizer(statusbar_sizer);

	m_bookCtrl = new wxNotebook(frame, wxID_ANY, wxDefaultPosition, wxSize(200, 100), 0);

	wxScrolled<wxPanel> *panelObject = new wxScrolled<wxPanel>(m_bookCtrl);	
	panelObject->EnableScrolling(true, true);
	panelObject->SetScrollbars(5, 10, 20, 60);
	wxBoxSizer * panelObject_sizer = new wxBoxSizer(wxVERTICAL);
	objectName = new wxTextCtrl(panelObject, 1000, "", wxDefaultPosition, wxDefaultSize);
	show_edges = new wxCheckBox(panelObject, EDGES_CHECKBOX, "Show Edges", wxDefaultPosition, wxDefaultSize);
	Connect(EDGES_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	interp_normals = new wxCheckBox(panelObject, INTERP_CHECKBOX, "Interp. Normals per Vtx", wxDefaultPosition, wxDefaultSize);
	Connect(INTERP_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);

	/*wxBoxSizer * transp_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* transp_text = new wxStaticText(panelObject, 9999, "Refraction index: ");
	transparent = new wxCheckBox(panelObject, TRANSPARENT_CHECKBOX, "Transparent", wxDefaultPosition, wxDefaultSize);
	Connect(TRANSPARENT_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	refractionIndex = new wxSpinCtrlDouble(panelObject, REFRACTION_INDEX, "1.5", wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0.1, 10, 1.5, 0.05);
	Connect(REFRACTION_INDEX, wxEVT_SPINCTRLDOUBLE, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	transp_sizer->Add(transparent);
	transp_sizer->Add(transp_text);
	transp_sizer->Add(refractionIndex);*/

	flipnormals = new wxCheckBox(panelObject, FLIPNORMALS_CHECKBOX, "Flip normals for transparency", wxDefaultPosition, wxDefaultSize);
	Connect(FLIPNORMALS_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);

	ghost = new wxCheckBox(panelObject, GHOST_CHECKBOX, "Ghost object", wxDefaultPosition, wxDefaultSize);
	Connect(GHOST_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);

	m_AlbedoFile = new wxListCtrl(panelObject, ALBEDO_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	wxListItem itCol;
	itCol.SetId(0);
	itCol.SetText("Albedo Texture");
	itCol.SetWidth(200);
	wxListItem itCol2;
	itCol2.SetId(1);
	itCol2.SetText("Albedo Color");
	itCol2.SetWidth(200);
	m_AlbedoFile->InsertColumn(0, itCol);
	m_AlbedoFile->InsertColumn(1, itCol2);
	m_AlbedoFile->SetDropTarget(new DnDAlbedoFile(m_AlbedoFile, frame));
	Connect(ALBEDO_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClick), NULL, frame);
	Connect(ALBEDO_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClick), NULL, frame);
	Connect(ALBEDO_FILES, wxEVT_LIST_ITEM_SELECTED, wxListEventHandler(RaytracerFrame::OnListSelected), NULL, frame);

	m_SpecularFile = new wxListCtrl(panelObject, SPECULAR_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	wxListItem itCol4;
	itCol4.SetId(0);
	itCol4.SetText("Specular Map");
	itCol4.SetWidth(200);
	itCol2.SetId(1);
	itCol2.SetText("Specular Color");
	itCol2.SetWidth(200);
	m_SpecularFile->InsertColumn(0, itCol4);
	m_SpecularFile->InsertColumn(1, itCol2);
	m_SpecularFile->SetDropTarget(new DnDSpecularFile(m_SpecularFile, frame));
	Connect(SPECULAR_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickSpecular), NULL, frame);
	Connect(SPECULAR_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickSpecular), NULL, frame);
	Connect(SPECULAR_FILES, wxEVT_LIST_ITEM_SELECTED, wxListEventHandler(RaytracerFrame::OnListSelected), NULL, frame);

	m_RoughnessFile = new wxListCtrl(panelObject, ROUGHNESS_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	wxListItem itCol5;
	itCol5.SetId(0);
	itCol5.SetText("Phong Exponent Map");
	itCol5.SetWidth(200);
	itCol2.SetId(1);
	itCol2.SetText("Exponent Values");
	itCol2.SetWidth(200);
	m_RoughnessFile->InsertColumn(0, itCol5);
	m_RoughnessFile->InsertColumn(1, itCol2);
	m_RoughnessFile->SetDropTarget(new DnDRoughnessFile(m_RoughnessFile, frame));
	Connect(ROUGHNESS_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickRoughness), NULL, frame);
	Connect(ROUGHNESS_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickRoughness), NULL, frame);
	Connect(ROUGHNESS_FILES, wxEVT_LIST_ITEM_SELECTED, wxListEventHandler(RaytracerFrame::OnListSelected), NULL, frame);

	m_NormalFile = new wxListCtrl(panelObject, NORMAL_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	itCol2.SetId(0);
	itCol2.SetText("Normal Map");
	itCol2.SetWidth(300);
	m_NormalFile->InsertColumn(0, itCol2);
	m_NormalFile->SetDropTarget(new DnDNormalFile(m_NormalFile, frame));
	Connect(NORMAL_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickNormal), NULL, frame);
	Connect(NORMAL_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickNormal), NULL, frame);
	Connect(NORMAL_FILES, wxEVT_LIST_ITEM_SELECTED, wxListEventHandler(RaytracerFrame::OnListSelected), NULL, frame);

	m_SubsurfaceFile = new wxListCtrl(panelObject, SUBSURFACE_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	wxListItem itCol6;
	itCol6.SetId(0);
	itCol6.SetText("Subsurface scattering Map");
	itCol6.SetWidth(200);
	itCol2.SetId(1);
	itCol2.SetText("Subsurface scattering Color");
	itCol2.SetWidth(200);
	m_SubsurfaceFile->InsertColumn(0, itCol6);
	m_SubsurfaceFile->InsertColumn(1, itCol2);
	m_SubsurfaceFile->SetDropTarget(new DnDSubsurfaceFile(m_SubsurfaceFile, frame));
	Connect(SUBSURFACE_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickSubsurface), NULL, frame);
	Connect(SUBSURFACE_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickSubsurface), NULL, frame);
	Connect(SUBSURFACE_FILES, wxEVT_LIST_ITEM_SELECTED, wxListEventHandler(RaytracerFrame::OnListSelected), NULL, frame);

	//albedoColorPicker = new wxColourPickerCtrl(panelObject, ALBEDO_COLORPICKER, wxColour(255,255,255), wxDefaultPosition, wxDefaultSize, wxCLRP_USE_TEXTCTRL | wxCLRP_SHOW_LABEL);
	//Connect(ALBEDO_COLORPICKER, wxEVT_COLOURPICKER_CHANGED, wxCommandEventHandler(RenderPanel<double>::update_parameters_and_render), NULL, renderPanel);

	m_AlphaFile = new wxListCtrl(panelObject, ALPHA_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	wxListItem itCol3;
	itCol3.SetId(0);
	itCol3.SetText("Alpha Map");
	itCol3.SetWidth(300);
	itCol2.SetId(1);
	itCol2.SetText("Transparent ?");
	itCol2.SetWidth(200);
	m_AlphaFile->InsertColumn(0, itCol3);
	m_AlphaFile->InsertColumn(1, itCol2);
	m_AlphaFile->SetDropTarget(new DnDAlphaFile(m_AlphaFile, frame));
	Connect(ALPHA_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickAlpha), NULL, frame);
	Connect(ALPHA_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickAlpha), NULL, frame);
	Connect(ALPHA_FILES, wxEVT_LIST_ITEM_SELECTED, wxListEventHandler(RaytracerFrame::OnListSelected), NULL, frame);


	m_TranspFile = new wxListCtrl(panelObject, TRANSP_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	wxListItem itCol7;
	itCol7.SetId(0);
	itCol7.SetText("Refractive Map");
	itCol7.SetWidth(200);
	itCol2.SetId(1);
	itCol2.SetText("Refractive ?");
	itCol2.SetWidth(200);
	m_TranspFile->InsertColumn(0, itCol7);
	m_TranspFile->InsertColumn(1, itCol2);
	m_TranspFile->SetDropTarget(new DnDTranspFile(m_TranspFile, frame));
	Connect(TRANSP_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickTransparent), NULL, frame);
	Connect(TRANSP_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickTransparent), NULL, frame);
	Connect(TRANSP_FILES, wxEVT_LIST_ITEM_SELECTED, wxListEventHandler(RaytracerFrame::OnListSelected), NULL, frame);


	m_RefrFile = new wxListCtrl(panelObject, REFR_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	wxListItem itCol8;
	itCol8.SetId(0);
	itCol8.SetText("Refraction Index Map");
	itCol8.SetWidth(200);
	itCol2.SetId(1);
	itCol2.SetText("Refraction Index");
	itCol2.SetWidth(200);
	m_RefrFile->InsertColumn(0, itCol8);
	m_RefrFile->InsertColumn(1, itCol2);
	m_RefrFile->SetDropTarget(new DnDRefrFile(m_RefrFile, frame));
	Connect(REFR_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickRefrIndex), NULL, frame);
	Connect(REFR_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickRefrIndex), NULL, frame);
	Connect(REFR_FILES, wxEVT_LIST_ITEM_SELECTED, wxListEventHandler(RaytracerFrame::OnListSelected), NULL, frame);

	/*wxBoxSizer * ks_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* ks_text = new wxStaticText(panelObject, 9999, "Ks : ");
	ks_slider = new wxSlider(panelObject, KS_SLIDER, 0, 0, 100, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL);
	Connect(KS_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel<double>::update_parameters_and_render), NULL, renderPanel);
	ks_sizer->Add(ks_text, 0, wxEXPAND);
	ks_sizer->Add(ks_slider, 1, wxEXPAND);*/

	colorAnisotropy = new wxButton(panelObject, COLOR_ANISOTROPY, "Color Anisotropy", wxDefaultPosition, wxDefaultSize);
	Connect(COLOR_ANISOTROPY, wxEVT_BUTTON, wxCommandEventHandler(RenderPanel::colorAnisotropy), NULL, renderPanel);

	randomColors = new wxButton(panelObject, RANDOM_COLOR, "Randomize Face Color Map", wxDefaultPosition, wxDefaultSize);
	Connect(RANDOM_COLOR, wxEVT_BUTTON, wxCommandEventHandler(RenderPanel::randomColors), NULL, renderPanel);

	deleteObject = new wxButton(panelObject, DELETE_OBJECT, "Delete", wxDefaultPosition, wxDefaultSize);
	Connect(DELETE_OBJECT, wxEVT_BUTTON, wxCommandEventHandler(RenderPanel::delete_object), NULL, renderPanel);

	panelObject_sizer->Add(objectName, 0, wxEXPAND);
	panelObject_sizer->Add(show_edges, 0, wxEXPAND);
	panelObject_sizer->Add(interp_normals, 0, wxEXPAND);
	//panelObject_sizer->Add(transp_sizer, 0, wxEXPAND);
	panelObject_sizer->Add(flipnormals, 0, wxEXPAND);
	panelObject_sizer->Add(ghost, 0, wxEXPAND);
	panelObject_sizer->Add(m_AlbedoFile, 0, wxEXPAND);
	//panelObject_sizer->Add(albedoColorPicker, 0, wxEXPAND);
	panelObject_sizer->Add(m_SpecularFile, 0, wxEXPAND);
	panelObject_sizer->Add(m_RoughnessFile, 0, wxEXPAND);
	panelObject_sizer->Add(m_NormalFile, 0, wxEXPAND);
	panelObject_sizer->Add(m_SubsurfaceFile, 0, wxEXPAND);
	panelObject_sizer->Add(m_AlphaFile, 0, wxEXPAND);
	panelObject_sizer->Add(m_TranspFile, 0, wxEXPAND);
	panelObject_sizer->Add(m_RefrFile, 0, wxEXPAND);
	//panelObject_sizer->Add(ks_sizer, 0, wxEXPAND);
	panelObject_sizer->Add(colorAnisotropy, 0, wxEXPAND);
	panelObject_sizer->Add(randomColors, 0, wxEXPAND);
	panelObject_sizer->Add(deleteObject, 0, wxEXPAND);
	
	panelObject->SetSizer(panelObject_sizer);


	m_bookCtrl->AddPage(panelObject, wxT("Object"), false);


	wxPanel *panelRenderer = new wxPanel(m_bookCtrl);
	wxBoxSizer * panelRenderer_sizer = new wxBoxSizer(wxVERTICAL);

	wxBoxSizer * rendersize_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* renderwidth_text = new wxStaticText(panelRenderer, 9999 - 1, "W: ");
#ifdef _DEBUG
	renderwidth = new wxSpinCtrl(panelRenderer, RENDER_WIDTH, wxString("100"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 2, 10000, 100);	
#else
	renderwidth = new wxSpinCtrl(panelRenderer, RENDER_WIDTH, wxString("1000"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 2, 10000, 1000);
#endif
	Connect(RENDER_WIDTH, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	wxStaticText* renderheight_text = new wxStaticText(panelRenderer, 9999 - 1, "H: ");
#ifdef _DEBUG
	renderheight = new wxSpinCtrl(panelRenderer, RENDER_HEIGHT, wxString("80"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 2, 10000, 80);
#else
	renderheight = new wxSpinCtrl(panelRenderer, RENDER_HEIGHT, wxString("800"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 2, 10000, 800);
#endif
	Connect(RENDER_HEIGHT, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	rendersize_sizer->Add(renderwidth_text, 0, wxEXPAND);
	rendersize_sizer->Add(renderwidth, 1, wxEXPAND);
	rendersize_sizer->Add(renderheight_text, 0, wxEXPAND);
	rendersize_sizer->Add(renderheight, 1, wxEXPAND);

	wxBoxSizer * nbrays_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* nbrays_text = new wxStaticText(panelRenderer, 9999 - 1, "nb paths per pixel: ");
	nbrays = new wxSpinCtrl(panelRenderer, NBRAYS, wxString("100"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 10000, 100);
	Connect(NBRAYS, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
#if USE_OPENIMAGEDENOISER
	has_denoiser = new wxCheckBox(panelRenderer, FILTER_CHECKBOX, "Denoiser (not temp. consistent)", wxDefaultPosition, wxDefaultSize);
	has_denoiser->SetToolTip("ONLY for offline rendering (no previz). Beware, not temporally stable.");
	Connect(FILTER_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
#endif
	nbrays_sizer->Add(nbrays_text, 0, wxEXPAND);
	nbrays_sizer->Add(nbrays, 1, wxEXPAND);
	nbrays_sizer->Add(has_denoiser, 1, wxEXPAND);

	wxBoxSizer * fov_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* fov_text = new wxStaticText(panelRenderer, 9999, "fov (deg): ");
	fov_slider = new wxSlider(panelRenderer, FOV_SLIDER, 35, 2, 179, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL|wxSL_LABELS);
	Connect(FOV_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	fov_sizer->Add(fov_text, 0, wxEXPAND);
	fov_sizer->Add(fov_slider, 1, wxEXPAND);

	wxBoxSizer * aperture_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* aperture_text = new wxStaticText(panelRenderer, 9999-1, "aperture (mm): ");
	aperture_slider = new wxSlider(panelRenderer, APERTURE_SLIDER, 10, 0, 3000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(APERTURE_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	aperture_sizer->Add(aperture_text, 0, wxEXPAND);
	aperture_sizer->Add(aperture_slider, 1, wxEXPAND);

	wxBoxSizer * focus_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* focus_text = new wxStaticText(panelRenderer, 9999 - 1, "focus distance (cm): ");
	focus_slider = new wxSlider(panelRenderer, FOCUS_SLIDER, 5000, 0, 10000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(FOCUS_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	focus_sizer->Add(focus_text, 0, wxEXPAND);
	focus_sizer->Add(focus_slider, 1, wxEXPAND);

	wxBoxSizer * filter_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* filter_text = new wxStaticText(panelRenderer, 9999 - 1, "filter width: ");
	filter_slider = new wxSlider(panelRenderer, FILTER_SLIDER, 5, 1, 100, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(FILTER_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	filter_sizer->Add(filter_text, 0, wxEXPAND);
	filter_sizer->Add(filter_slider, 1, wxEXPAND);


	wxBoxSizer * envmap_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* envmap_text = new wxStaticText(panelRenderer, 9999 - 1, "envmap: ");
	envmapName = new wxTextCtrl(panelRenderer, 1000, "", wxDefaultPosition, wxDefaultSize);
	envmapName->SetDropTarget(new DnDEnvmapFile(envmapName, frame));
	envmap_sizer->Add(envmap_text, 0, wxEXPAND);
	envmap_sizer->Add(envmapName, 1, wxEXPAND);

	wxBoxSizer * envmapintensity_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* envmapintensity_text = new wxStaticText(panelRenderer, 9999 - 1, "envmap intensity: ");
	envmapintensity_slider = new wxSlider(panelRenderer, ENVMAPINTENSITY_SLIDER, 100, 0, 10000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(ENVMAPINTENSITY_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	envmapintensity_sizer->Add(envmapintensity_text, 0, wxEXPAND);
	envmapintensity_sizer->Add(envmapintensity_slider, 1, wxEXPAND);

	wxBoxSizer * background_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* background_text = new wxStaticText(panelRenderer, 9999 - 1, "backgrd image: ");
	backgroundName = new wxTextCtrl(panelRenderer, 1000, "", wxDefaultPosition, wxDefaultSize);
	backgroundName->SetDropTarget(new DnDBackgroundFile(backgroundName, frame));
	background_sizer->Add(background_text, 0, wxEXPAND);
	background_sizer->Add(backgroundName, 1, wxEXPAND);

	wxBoxSizer * lightintensity_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* lightintensity_text = new wxStaticText(panelRenderer, 9999 - 1, "light intensity: ");
	lightintensity_slider = new wxSlider(panelRenderer, LIGHTINTENSITY_SLIDER, 100, 0, 1000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(LIGHTINTENSITY_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	lightintensity_sizer->Add(lightintensity_text, 0, wxEXPAND);
	lightintensity_sizer->Add(lightintensity_slider, 1, wxEXPAND);

	wxBoxSizer * bounces_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* bounces_text = new wxStaticText(panelRenderer, 9999 - 1, "light bounces: ");
	bounces = new wxSpinCtrl(panelRenderer, BOUNCES_SPIN, wxString("3"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 100, 3);
	Connect(BOUNCES_SPIN, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	bounces_sizer->Add(bounces_text, 0, wxEXPAND);
	bounces_sizer->Add(bounces, 1, wxEXPAND);

	launchRender = new wxButton(panelRenderer, LAUNCH_RENDER, "Offline render (freezes app)", wxDefaultPosition, wxDefaultSize);
	Connect(LAUNCH_RENDER, wxEVT_BUTTON, wxCommandEventHandler(RenderPanel::launch_render), NULL, renderPanel);
	
	panelRenderer_sizer->Add(rendersize_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(nbrays_sizer, 0, wxEXPAND);	
	panelRenderer_sizer->Add(fov_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(aperture_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(focus_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(filter_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(envmap_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(envmapintensity_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(background_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(lightintensity_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(bounces_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(launchRender, 0, wxEXPAND);
	panelRenderer->SetSizer(panelRenderer_sizer);

	m_bookCtrl->AddPage(panelRenderer, wxT("Renderer"), false);



	wxPanel *panelFog = new wxPanel(m_bookCtrl);
	wxBoxSizer * panelFog_sizer = new wxBoxSizer(wxVERTICAL);

	wxBoxSizer * fogdensity_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* fogdensity_text = new wxStaticText(panelFog, 9999 - 1, "fog density : ");
	fogdensity_slider = new wxSlider(panelFog, FOGDENSITY_SLIDER, 0, 0, 100, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(FOGDENSITY_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	fogdensity_sizer->Add(fogdensity_text, 0, wxEXPAND);
	fogdensity_sizer->Add(fogdensity_slider, 1, wxEXPAND);

	uniformFogRadio = new wxRadioButton(panelFog, UNIFORMFOG_RADIO, "Uniform fog", wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
	expFogRadio = new wxRadioButton(panelFog, EXPFOG_RADIO, "Exponential fog", wxDefaultPosition, wxDefaultSize, 0);
	Connect(UNIFORMFOG_RADIO, wxEVT_RADIOBUTTON, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	Connect(UNIFORMFOG_RADIO, wxEVT_RADIOBUTTON, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);

	panelFog_sizer->Add(fogdensity_sizer, 0, wxEXPAND);
	panelFog_sizer->Add(uniformFogRadio, 0, wxEXPAND);
	panelFog_sizer->Add(expFogRadio, 0, wxEXPAND);
	panelFog->SetSizer(panelFog_sizer);

	m_bookCtrl->AddPage(panelFog, wxT("Fog"), false);



	wxPanel *panelLenticular = new wxPanel(m_bookCtrl);
	wxBoxSizer * panelLenticular_sizer = new wxBoxSizer(wxVERTICAL);
	isLenticularCheck = new wxCheckBox(panelLenticular, IS_LENTICULAR_CHECKBOX, "Lenticular Camera", wxDefaultPosition, wxDefaultSize);
	Connect(IS_LENTICULAR_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	panelLenticular_sizer->Add(isLenticularCheck, 0, wxEXPAND);

	wxBoxSizer * nbviews_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* nbviews_text = new wxStaticText(panelLenticular, 9999 - 1, "nb views: ");
	nbviews = new wxSpinCtrl(panelLenticular, NBVIEWS, wxString("10"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 200, 10);
	Connect(NBVIEWS, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	nbviews_sizer->Add(nbviews_text, 0, wxEXPAND);
	nbviews_sizer->Add(nbviews, 1, wxEXPAND);
	panelLenticular_sizer->Add(nbviews_sizer, 0, wxEXPAND);

	wxBoxSizer * nbpixslice_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* nbpixslice_text = new wxStaticText(panelLenticular, 9999 - 1, "nb pixels per slice: ");
	nbpixslice = new wxSpinCtrl(panelLenticular, NBPIXPERSLICE, wxString("1"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 200, 1);
	Connect(NBPIXPERSLICE, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	nbpixslice_sizer->Add(nbpixslice_text, 0, wxEXPAND);
	nbpixslice_sizer->Add(nbpixslice, 1, wxEXPAND);
	panelLenticular_sizer->Add(nbpixslice_sizer, 0, wxEXPAND);

	wxBoxSizer * maxangle_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* maxangle_text = new wxStaticText(panelLenticular, 9999 - 1, "max disparity (angle) : ");
	maxangle_slider = new wxSlider(panelLenticular, MAXANGLE_SLIDER, 38, 0, 180, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(MAXANGLE_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	maxangle_sizer->Add(maxangle_text, 0, wxEXPAND);
	maxangle_sizer->Add(maxangle_slider, 1, wxEXPAND);
	panelLenticular_sizer->Add(maxangle_sizer, 0, wxEXPAND);


	isArrayCheck = new wxCheckBox(panelLenticular, IS_ARRAY_CHECKBOX, "Camera Array", wxDefaultPosition, wxDefaultSize);
	Connect(IS_ARRAY_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	panelLenticular_sizer->Add(isArrayCheck, 0, wxEXPAND);

	wxBoxSizer * nbviewsXY_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* nbviewsXY_text = new wxStaticText(panelLenticular, 9999 - 1, "nb views x, y: ");
	nbviewsX = new wxSpinCtrl(panelLenticular, NBVIEWSX, wxString("1"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 1000, 1);
	Connect(NBVIEWSX, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	nbviewsY = new wxSpinCtrl(panelLenticular, NBVIEWSY, wxString("1"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 1000, 1);
	Connect(NBVIEWSY, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	nbviewsXY_sizer->Add(nbviewsXY_text, 0, wxEXPAND);
	nbviewsXY_sizer->Add(nbviewsX, 1, wxEXPAND);
	nbviewsXY_sizer->Add(nbviewsY, 1, wxEXPAND);
	panelLenticular_sizer->Add(nbviewsXY_sizer, 0, wxEXPAND);

	wxBoxSizer * spacing_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* spacing_text = new wxStaticText(panelLenticular, 9999 - 1, "Max distance x, y: ");
	maxspacingX = new wxSpinCtrlDouble(panelLenticular, MAXSPACINGX, wxString("3.0"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0.01, 200, 1.0, 0.1);
	Connect(MAXSPACINGX, wxEVT_SPINCTRLDOUBLE, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	maxspacingY = new wxSpinCtrlDouble(panelLenticular, MAXSPACINGY, wxString("3.0"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0.01, 200, 1.0, 0.1);
	Connect(MAXSPACINGY, wxEVT_SPINCTRLDOUBLE, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	spacing_sizer->Add(spacing_text, 0, wxEXPAND);
	spacing_sizer->Add(maxspacingX, 1, wxEXPAND);
	spacing_sizer->Add(maxspacingY, 1, wxEXPAND);
	panelLenticular_sizer->Add(spacing_sizer, 0, wxEXPAND);


	wxBoxSizer * doublefrustum_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* doublefrustum_text = new wxStaticText(panelLenticular, 9999 - 1, "double frustum (ray start): ");
	doubleFrustum = new wxSpinCtrlDouble(panelLenticular, DOUBLEFRUSTUM, wxString("0.0"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, -200, 200, 0.0, 0.1);
	Connect(DOUBLEFRUSTUM, wxEVT_SPINCTRLDOUBLE, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	doublefrustum_sizer->Add(doublefrustum_text, 0, wxEXPAND);
	doublefrustum_sizer->Add(doubleFrustum, 1, wxEXPAND);
	panelLenticular_sizer->Add(doublefrustum_sizer, 0, wxEXPAND);
	

	panelLenticular->SetSizer(panelLenticular_sizer);
	m_bookCtrl->AddPage(panelLenticular, wxT("Lenticular/Holography"), false);






	wxPanel *panelAnimation = new wxPanel(m_bookCtrl);
	wxBoxSizer * panelAnimation_sizer = new wxBoxSizer(wxVERTICAL);
	panelAnimation->SetSizer(panelAnimation_sizer);

	wxBoxSizer * nbframes_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* nbframes_text = new wxStaticText(panelAnimation, 9999 - 1, "nb frames: ");
	nbframesctrl = new wxSpinCtrl(panelAnimation, NBFRAMES, wxString("1"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 10000, 1);
	Connect(NBFRAMES, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	nbframes_sizer->Add(nbframes_text, 0, wxEXPAND);
	nbframes_sizer->Add(nbframesctrl, 1, wxEXPAND);

	wxBoxSizer * duration_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* duration_text = new wxStaticText(panelAnimation, 9999 - 1, "duration (s): ");
	duration = new wxSpinCtrlDouble(panelAnimation, DURATION, wxString("1"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0.001, 10000, 1, 0.1);
	Connect(DURATION, wxEVT_SPINCTRLDOUBLE, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	duration_sizer->Add(duration_text, 0, wxEXPAND);
	duration_sizer->Add(duration, 1, wxEXPAND);

	panelAnimation_sizer->Add(nbframes_sizer, 0, wxEXPAND);
	panelAnimation_sizer->Add(duration_sizer, 0, wxEXPAND);

	recordKeyframes = new wxToggleButton(panelAnimation, RECORD_KEYFRAME, "Enable Record keyframes", wxDefaultPosition, wxDefaultSize);
	Connect(RECORD_KEYFRAME, wxEVT_TOGGLEBUTTON, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	panelAnimation_sizer->Add(recordKeyframes, 0, wxEXPAND);

	addKeyframe = new wxButton(panelAnimation, ADD_KEYFRAME, "Add object keyframe", wxDefaultPosition, wxDefaultSize);
	Connect(ADD_KEYFRAME, wxEVT_BUTTON, wxCommandEventHandler(RenderPanel::add_keyframe), NULL, renderPanel);
	panelAnimation_sizer->Add(addKeyframe, 0, wxEXPAND);

	renderVideo = new wxButton(panelAnimation, RENDER_VIDEO, "Offline render video (freezes)", wxDefaultPosition, wxDefaultSize);
	Connect(RENDER_VIDEO, wxEVT_BUTTON, wxCommandEventHandler(RenderPanel::render_video), NULL, renderPanel);
	panelAnimation_sizer->Add(renderVideo, 0, wxEXPAND);

	m_bookCtrl->AddPage(panelAnimation, wxT("Animation"), false);

	// ------------------
	wxPanel *panelFluids = new wxPanel(m_bookCtrl);
	wxBoxSizer * panelFluids_sizer = new wxBoxSizer(wxVERTICAL);
	panelFluids->SetSizer(panelFluids_sizer);

	wxBoxSizer * fluidres_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* fluidres_text = new wxStaticText(panelFluids, 9999 - 1, "Resolution: ");
	fluidresX = new wxSpinCtrl(panelFluids, FLUIDRESX, wxString("64"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 4, 512, 32);
	fluidresY = new wxSpinCtrl(panelFluids, FLUIDRESY, wxString("64"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 4, 512, 32);
	fluidresZ = new wxSpinCtrl(panelFluids, FLUIDRESZ, wxString("64"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 4, 512, 32);
	Connect(FLUIDRESX, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	Connect(FLUIDRESY, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	Connect(FLUIDRESZ, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	fluidres_sizer->Add(fluidres_text, 0, wxEXPAND);
	fluidres_sizer->Add(fluidresX, 1, wxEXPAND);
	fluidres_sizer->Add(fluidresY, 1, wxEXPAND);
	fluidres_sizer->Add(fluidresZ, 1, wxEXPAND);
	panelFluids_sizer->Add(fluidres_sizer, 0, wxEXPAND);

	wxBoxSizer * fluidnparticles_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* fluidnparticles_text = new wxStaticText(panelFluids, 9999 - 1, "Nb particles: ");
	fluidnparticles = new wxSpinCtrl(panelFluids, FLUIDNPARTICLES, wxString("1000000"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 10000000, 4000000);
	Connect(FLUIDNPARTICLES, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	fluidnparticles_sizer->Add(fluidnparticles_text, 0, wxEXPAND);
	fluidnparticles_sizer->Add(fluidnparticles, 1, wxEXPAND);
	panelFluids_sizer->Add(fluidnparticles_sizer, 0, wxEXPAND);

	wxBoxSizer * fluidparticlesize_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* fluidparticlesize_text = new wxStaticText(panelFluids, 9999 - 1, "Particle size: ");
	fluidparticlesize = new wxSpinCtrlDouble(panelFluids, FLUIDPARTICLESIZE, wxString("0.16"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0.01, 1, 0.05, 0.01);
	Connect(FLUIDPARTICLESIZE, wxEVT_SPINCTRLDOUBLE, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	fluidparticlesize_sizer->Add(fluidparticlesize_text, 0, wxEXPAND);
	fluidparticlesize_sizer->Add(fluidparticlesize, 1, wxEXPAND);
	panelFluids_sizer->Add(fluidparticlesize_sizer, 0, wxEXPAND);

	wxBoxSizer * fluidtimestep_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* fluidtimestep_text = new wxStaticText(panelFluids, 9999 - 1, "Time step: ");
	fluidtimestep = new wxSpinCtrlDouble(panelFluids, FLUIDTIMESTEP, wxString("0.025"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0.001, 2, 1./30., 0.01);
	Connect(FLUIDTIMESTEP, wxEVT_SPINCTRLDOUBLE, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	fluidtimestep_sizer->Add(fluidtimestep_text, 0, wxEXPAND);
	fluidtimestep_sizer->Add(fluidtimestep, 1, wxEXPAND);
	panelFluids_sizer->Add(fluidtimestep_sizer, 0, wxEXPAND);


	wxBoxSizer * fluidsubsteps_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* fluidsubsteps_text = new wxStaticText(panelFluids, 9999 - 1, "Nb substeps: ");
	fluidsubsteps = new wxSpinCtrl(panelFluids, FLUIDNSUBSTEPS, wxString("5"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 1000, 5);
	Connect(FLUIDNSUBSTEPS, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	fluidsubsteps_sizer->Add(fluidsubsteps_text, 0, wxEXPAND);
	fluidsubsteps_sizer->Add(fluidsubsteps, 1, wxEXPAND);
	panelFluids_sizer->Add(fluidsubsteps_sizer, 0, wxEXPAND);
	

	wxBoxSizer *initfluid_sizer = new wxBoxSizer(wxHORIZONTAL);	
	initfluid = new wxCheckBox(panelFluids, IS_FLUID_INIT_CHECKBOX, "Init fluid with selected mesh", wxDefaultPosition, wxDefaultSize);
	Connect(IS_FLUID_INIT_CHECKBOX, wxEVT_CHECKBOX, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);	
	initfluid_sizer->Add(initfluid, 1, wxEXPAND);
	panelFluids_sizer->Add(initfluid_sizer, 0, wxEXPAND);

	

	addFluid = new wxButton(panelFluids, ADD_FLUID, "Add fluid", wxDefaultPosition, wxDefaultSize);
	Connect(ADD_FLUID, wxEVT_BUTTON, wxCommandEventHandler(RenderPanel::add_fluid), NULL, renderPanel);
	panelFluids_sizer->Add(addFluid, 0, wxEXPAND);



	
	m_bookCtrl->AddPage(panelFluids, wxT("Fluids"), false);

	// ----------------

	m_bookCtrl->SetSelection(0);

	wxBoxSizer * renderSuperSizer_sizer = new wxBoxSizer(wxVERTICAL);
	renderSuperSizer_sizer->Add(renderPanel,5, wxEXPAND);

	wxBoxSizer * time_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* time_text = new wxStaticText(frame, 9999, "timeline ");
	time_slider = new wxSlider(frame, FOV_SLIDER, 0, 0, 1, wxDefaultPosition, wxSize(600,32), wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(TIME_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);	
	time_sizer->Add(time_text, 0, wxEXPAND);
	time_sizer->Add(time_slider, 1, wxEXPAND);

	renderSuperSizer_sizer->Add(time_sizer);

	wxBoxSizer * sizer = new wxBoxSizer(wxHORIZONTAL);
	frame->SetSizer(sizer);
	sizer->Add(renderSuperSizer_sizer, 2, wxEXPAND);
	sizer->Add(m_bookCtrl, 1, wxEXPAND);




	colPicker = new wxColourDialog(frame);
	texOpenDlg = new wxFileDialog(frame, _("Open texture file"), "", "", "All files (*.*)|*.*", wxFD_OPEN | wxFD_FILE_MUST_EXIST);

	render_loop_on = false;
	frame->Show();
	activateRenderLoop(true);

    return true;
}

void RaytracerFrame::OnSize(wxSizeEvent& event) {
	if (render_panel) {
		int displayWidth = render_panel->GetClientSize().GetX();
		render_panel->displayH = std::min(render_panel->GetClientSize().GetY()-1, (int)(render_panel->raytracer.H / (double)render_panel->raytracer.W * displayWidth));
		render_panel->displayW = std::min(render_panel->GetClientSize().GetX() - 1, (int)(render_panel->displayH * render_panel->raytracer.W/ (double)render_panel->raytracer.H));
		render_panel->displayH = std::min(render_panel->GetClientSize().GetY() - 1, (int)(render_panel->raytracer.H / (double)render_panel->raytracer.W * render_panel->displayW));
	}
	Refresh();

	event.Skip();
}

void RaytracerApp::activateRenderLoop(bool on)
{
#if !defined(__WXOSX__)&&!defined(__WXMAC__)&&!defined(__APPLE__)
	renderPanel->SetDoubleBuffered(true);
#endif
	if (on && !render_loop_on)
	{
		Connect(wxID_ANY, wxEVT_IDLE, wxIdleEventHandler(RaytracerApp::OnIdle));
		Connect(wxID_ANY, wxEVT_ERASE_BACKGROUND, wxEraseEventHandler(RaytracerApp::OnEraseBackGround));
		render_loop_on = true;
	} else if (!on && render_loop_on) {
		Disconnect(wxEVT_IDLE, wxIdleEventHandler(RaytracerApp::OnIdle));
		render_loop_on = false;
	}
}

RaytracerFrame::RaytracerFrame()
        : wxFrame(NULL, wxID_ANY, wxT("Pathtracer by N.Bonneel"),
                  wxPoint(10, 100))

{
    // frame icon and status bar
    //SetIcon(wxICON(sample));

    CreateStatusBar();
	programHandling = false;

	createAlbedoMenu();

    // construct menu
    wxMenu *file_menu = new wxMenu();
    file_menu->Append(ID_SAVE, wxT("&Save scene as..."));
	Connect(ID_SAVE, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::SaveAs), NULL, this);
	file_menu->Append(ID_OPEN, wxT("&Open scene..."));
	Connect(ID_OPEN, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::Open), NULL, this);
	file_menu->Append(ID_SAVEIMAGE, wxT("&Save current image..."));
	Connect(ID_SAVEIMAGE, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::SaveImage), NULL, this);
	file_menu->Append(ID_EXPORT_MTL, wxT("&Export .mtl (selected object)..."));
	Connect(ID_EXPORT_MTL, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::ExportMtl), NULL, this);

	wxMenu *info_menu = new wxMenu();
	info_menu->Append(ID_MESHINFO, wxT("&Mesh information..."));
	Connect(ID_MESHINFO, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::ShowMeshInfo), NULL, this);

	wxMenuBar* menu_bar = new wxMenuBar();
	menu_bar->Append(file_menu, wxT("&File"));	
	menu_bar->Append(info_menu, wxT("&Info"));

    SetMenuBar(menu_bar);

    SetBackgroundColour(*wxWHITE); // labels read better on this background
	
    this->SetDropTarget(new DnDFile(this));	

	this->SetSize(1024, 1024);
	Centre();

}




void RenderPanel::update_parameters_and_render(wxCommandEvent& event) {

	stop_render();

	raytracer.W = raytracer_app->renderwidth->GetValue();
	raytracer.H = raytracer_app->renderheight->GetValue();

	int displayWidth = GetClientSize().GetX();
	displayH = std::min(GetClientSize().GetY() - 1, (int)(raytracer.H / (double)raytracer.W * displayWidth));
	displayW = std::min(GetClientSize().GetX() - 1, (int)(displayH * raytracer.W / (double)raytracer.H));
	displayH = std::min(GetClientSize().GetY() - 1, (int)(raytracer.H / (double)raytracer.W * displayW));

	raytracer.nrays = raytracer_app->nbrays->GetValue();
	raytracer.cam.fov = raytracer_app->fov_slider->GetValue() * M_PI / 180.;
	raytracer.cam.aperture = raytracer_app->aperture_slider->GetValue() / 1000.;
	raytracer.cam.focus_distance = raytracer_app->focus_slider->GetValue() / 100.;
	raytracer.cam.is_lenticular = raytracer_app->isLenticularCheck->IsChecked();
	raytracer.sigma_filter = raytracer_app->filter_slider->GetValue() / 10.;

	raytracer.s.fog_density = raytracer_app->fogdensity_slider->GetValue() / 100.;

	raytracer.s.intensite_lumiere = raytracer_app->lightintensity_slider->GetValue() / 100. * 1000000000 * 4.*M_PI / (4.*M_PI*raytracer.s.lumiere->R*raytracer.s.lumiere->R*M_PI);
	raytracer.s.envmap_intensity = raytracer_app->envmapintensity_slider->GetValue() / 100.;
	raytracer.cam.lenticular_max_angle = raytracer_app->maxangle_slider->GetValue()*M_PI / 180;
	raytracer.cam.lenticular_nb_images = raytracer_app->nbviews->GetValue();
	raytracer.cam.lenticular_pixel_width = raytracer_app->nbpixslice->GetValue();

	raytracer.s.fog_type = raytracer_app->uniformFogRadio->GetValue() ? 0 : 1;
	raytracer.nb_bounces = raytracer_app->bounces->GetValue();

	raytracer.s.nbframes = raytracer_app->nbframesctrl->GetValue();
	raytracer_app->time_slider->SetMax(raytracer.s.nbframes);
	raytracer.s.duration = raytracer_app->duration->GetValue();
	raytracer.s.current_frame = raytracer_app->time_slider->GetValue();
	raytracer.s.current_time = raytracer.s.current_frame / (double)raytracer.s.nbframes*raytracer.s.duration;
	raytracer.is_recording = raytracer_app->recordKeyframes->GetValue();
	raytracer.s.double_frustum_start_t = raytracer_app->doubleFrustum->GetValue();

	raytracer.has_denoiser = raytracer_app->has_denoiser->GetValue();

	if (selected_object < 0 || (selected_object >= raytracer.s.objects.size())) {
		start_render();
		return;
	}


	raytracer.s.objects[selected_object]->display_edges = raytracer_app->show_edges->IsChecked();
	raytracer.s.objects[selected_object]->interp_normals = raytracer_app->interp_normals->IsChecked();	
	
	//wxColour alb = raytracer_app->albedoColorPicker->GetColour();
	//raytracer.s.objects[selected_object]->albedo = Vector(alb.Red()/255., alb.Green()/255., alb.Blue()/255.);

	//raytracer.s.objects[selected_object]->ks = raytracer_app->ks_slider->GetValue() / 100.;

	raytracer.s.objects[selected_object]->flip_normals = raytracer_app->flipnormals->IsChecked();	
	raytracer.s.objects[selected_object]->ghost = raytracer_app->ghost->IsChecked();
#ifdef USE_EMBREE
	if (raytracer.s.objects[selected_object]->type == OT_TRIMESH) {
		TriMesh* g = dynamic_cast<TriMesh*>(raytracer.s.objects[selected_object]);
		if (raytracer.s.objects[selected_object]->ghost)
			rtcSetGeometryMask(g->instance_geom, 1);
		else
			rtcSetGeometryMask(g->instance_geom, -1);
	}
#endif



	for (int i = 0; i < raytracer.s.objects.size(); i++) {
		raytracer.s.objects[i]->max_translation = raytracer.s.objects[i]->get_translation(raytracer.s.current_frame, false);
		raytracer.s.objects[i]->mat_rotation = raytracer.s.objects[i]->get_rotation(raytracer.s.current_frame, false);
		raytracer.s.objects[i]->scale = raytracer.s.objects[i]->get_scale(raytracer.s.current_frame, false);
	}
	for (int i = 0; i < raytracer.s.objects.size(); i++) {
		if (raytracer.s.objects[i]->type != OT_FLUID) continue;
		Fluid* f = dynamic_cast<Fluid*>(raytracer.s.objects[i]);
		f->build_bvh(raytracer_app->time_slider->GetValue(), 0, f->Nparticles);
		f->build_grid(raytracer_app->time_slider->GetValue());
	}

	if (selected_object >=0 && selected_object< raytracer.s.objects.size() && (raytracer.s.objects[selected_object]->type == OT_FLUID)) {
		dynamic_cast<Fluid*>(raytracer.s.objects[selected_object])->radius = raytracer_app->fluidparticlesize->GetValue();
	}


	start_render();
}
void RenderPanel::update_textures_and_render(wxCommandEvent& event) {

	if (selected_object < 0) return;
	if (selected_object >= raytracer.s.objects.size()) return;

	stop_render();

	raytracer.s.objects[selected_object]->display_edges = raytracer_app->show_edges->IsChecked();

	start_render();
}


void RenderPanel::add_fluid(wxCommandEvent& event) {
	raytracer_app->renderPanel->stop_render();

	Fluid* fl = new Fluid(raytracer.s, BBox(Vector(-18, -27.3, -18), Vector(18, -27.3+50, 18)), raytracer_app->fluidresX->GetValue(), raytracer_app->fluidresY->GetValue(), raytracer_app->fluidresZ->GetValue(), raytracer_app->fluidnparticles->GetValue(), 1000., raytracer_app->fluidparticlesize->GetValue(), raytracer_app->nbframesctrl->GetValue(), raytracer_app->fluidtimestep->GetValue(), raytracer_app->fluidsubsteps->GetValue());
	fl->init_particles(this->raytracer_app->initfluid->GetValue(), this->raytracer_app->renderPanel->selected_object);
	fl->run();
	
	raytracer.s.addObject(fl);
	raytracer_app->renderPanel->update_gui();
	raytracer_app->renderPanel->start_render();	
}

void RenderPanel::render_video(wxCommandEvent& event) {
	stop_render();
	PerfChrono perf;
	perf.Start();
	raytracer.stopped = false;

	for (int i = 0; i < raytracer.s.nbframes; i++) {
		raytracer.clear_image();
		raytracer.s.current_frame = i;
		raytracer.s.current_time = i * raytracer.s.duration / (double)raytracer.s.nbframes;

		for (int j = 0; j < raytracer.s.objects.size(); j++) {
			Fluid* f = dynamic_cast<Fluid*>(raytracer.s.objects[j]);
			if (!f) continue;
			f->build_bvh(i, 0, f->Nparticles);
			f->build_grid(i);
		}

		if (raytracer_app->isArrayCheck->GetValue()) {
			raytracer.cam.nbviewX = raytracer_app->nbviewsX->GetValue();
			raytracer.cam.nbviewY = raytracer_app->nbviewsY->GetValue();
			Vector pos = raytracer.cam.position;
			Vector right = cross(raytracer.cam.direction, raytracer.cam.up);
			double dx = raytracer_app->maxspacingX->GetValue() / raytracer_app->nbviewsX->GetValue();
			double dy = raytracer_app->maxspacingY->GetValue() / raytracer_app->nbviewsY->GetValue();
			raytracer.cam.isArray = true;
			for (int j = 0; j < raytracer_app->nbviewsY->GetValue(); j++) {
				for (int k = 0; k < raytracer_app->nbviewsX->GetValue(); k++) {
					raytracer.cam.current_viewX = k;
					raytracer.cam.current_viewY = j;

					raytracer.cam.position = pos + (k - raytracer.cam.nbviewX / 2)*dx * right + (-j+ raytracer.cam.nbviewY / 2)*dy * raytracer.cam.up;
					raytracer.clear_image();
					raytracer.render_image_nopreviz();
				}
			}
			raytracer.cam.position = pos;
		}
		else {
			raytracer.cam.isArray = false;
			raytracer.render_image_nopreviz();
		}

		
	}
	raytracer.stopped = true;
	start_render();
}

void RenderPanel::update_gui() {

	raytracer_app->renderwidth->SetValue(raytracer.W);
	raytracer_app->renderheight->SetValue(raytracer.H);

	raytracer_app->nbrays->SetValue(raytracer.nrays);
	raytracer_app->isLenticularCheck->SetValue(raytracer.cam.is_lenticular);
	raytracer_app->maxangle_slider->SetValue(raytracer.cam.lenticular_max_angle*180. / M_PI);
	raytracer_app->fov_slider->SetValue(raytracer.cam.fov * 180. / M_PI);
	raytracer_app->aperture_slider->SetValue(raytracer.cam.aperture * 1000);
	raytracer_app->focus_slider->SetValue(raytracer.cam.focus_distance*100.);
	raytracer_app->filter_slider->SetValue(raytracer.sigma_filter*10.);
	raytracer_app->bounces->SetValue(raytracer.nb_bounces);
	raytracer_app->nbviews->SetValue(raytracer.cam.lenticular_nb_images);
	raytracer_app->nbpixslice->SetValue(raytracer.cam.lenticular_pixel_width);
	raytracer_app->nbframesctrl->SetValue(raytracer.s.nbframes);
	raytracer_app->time_slider->SetMax(raytracer.s.nbframes);


	raytracer_app->fogdensity_slider->SetValue(raytracer.s.fog_density*100.);
	raytracer_app->uniformFogRadio->SetValue(raytracer.s.fog_type == 0);

	double factor = 1. / 100. * 1000000000 * 4.*M_PI / (4.*M_PI*raytracer.s.lumiere->R*raytracer.s.lumiere->R*M_PI);
	raytracer_app->lightintensity_slider->SetValue(raytracer.s.intensite_lumiere / factor);
	raytracer_app->envmapintensity_slider->SetValue(raytracer.s.envmap_intensity * 100);

	raytracer_app->duration->SetValue(raytracer.s.duration);
	raytracer_app->time_slider->SetValue(raytracer.s.current_frame);
	raytracer_app->recordKeyframes->SetValue(raytracer.is_recording);


	if ((selected_object < 0) || (selected_object >= raytracer.s.objects.size())) {
		raytracer_app->objectName->SetLabelText(" ");
		return;
	}

	
	if (raytracer.s.objects[selected_object]->type == OT_TRIMESH) {
		TriMesh* g = dynamic_cast<TriMesh*>(raytracer.s.objects[selected_object]);
		std::ostringstream os;
		os << "triangles: " << g->indices.size() << ", vertices: " << g->vertices.size() << ", bvh leaves size: " << g->max_bvh_triangles << ", bvh max depth: " << g->bvh_depth << ", bvh avg depth:" << g->bvh_avg_depth << ", bvh nb nodes: " << g->bvh_nb_nodes << ", mouse distance: " << selected_object_t;
		raytracer_app->infoModel->SetLabelText(os.str().c_str());
	} else {
		std::ostringstream os;
		os << "mouse distance: " << selected_object_t;
		raytracer_app->infoModel->SetLabelText(os.str().c_str());
	}


	raytracer_app->objectName->SetLabelText(raytracer.s.objects[selected_object]->name);
	raytracer_app->show_edges->SetValue(raytracer.s.objects[selected_object]->display_edges);
	raytracer_app->interp_normals->SetValue(raytracer.s.objects[selected_object]->interp_normals);
	//raytracer_app->ks_slider->SetValue(raytracer.s.objects[selected_object]->ks*100.);

	raytracer_app->flipnormals->SetValue(raytracer.s.objects[selected_object]->flip_normals);	
	raytracer_app->ghost->SetValue(raytracer.s.objects[selected_object]->ghost);


	//Vector alb = raytracer.s.objects[selected_object]->albedo;
	//raytracer_app->albedoColorPicker->SetColour(wxColour(alb[0] * 255., alb[1] * 255., alb[2] * 255.));

	raytracer_app->m_AlbedoFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->textures.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->textures[i].filename);
		wxListItem item;
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->textures[i].filename;
		long index = raytracer_app->m_AlbedoFile->InsertItem(i, item);
		raytracer_app->m_AlbedoFile->SetItem(index, 0, txt, -1);
		raytracer_app->m_AlbedoFile->SetItem(index, 1, raytracer.s.objects[selected_object]->textures[i].multiplier.toColorStr(), -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_AlbedoFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}

	raytracer_app->m_NormalFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->normal_map.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->normal_map[i].filename);
		wxListItem item;
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->normal_map[i].filename;
		long index = raytracer_app->m_NormalFile->InsertItem(i, item);
		raytracer_app->m_NormalFile->SetItem(index, 0, txt, -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_NormalFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}


	raytracer_app->m_SubsurfaceFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->subsurface.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->subsurface[i].filename);
		wxListItem item;
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->subsurface[i].filename;
		long index = raytracer_app->m_SubsurfaceFile->InsertItem(i, item);
		raytracer_app->m_SubsurfaceFile->SetItem(index, 0, txt, -1);		
		raytracer_app->m_SubsurfaceFile->SetItem(index, 1, raytracer.s.objects[selected_object]->subsurface[i].multiplier.toColorStr(), -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_SubsurfaceFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}

	if (raytracer.s.objects.size()>1 && dynamic_cast<Sphere*>(raytracer.s.objects[1]))
		raytracer_app->envmapName->SetLabelText(((Sphere*)raytracer.s.objects[1])->envmapfilename);
	raytracer_app->backgroundName->SetLabelText(raytracer.s.backgroundfilename);

	raytracer_app->m_AlphaFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->alphamap.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->alphamap[i].filename);
		wxListItem item;
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->alphamap[i].filename;
		long index = raytracer_app->m_AlphaFile->InsertItem(i, item);
		raytracer_app->m_AlphaFile->SetItem(index, 0, txt, -1);
		raytracer_app->m_AlphaFile->SetItem(index, 1, (raytracer.s.objects[selected_object]->alphamap[i].multiplier[0]<0.5) ? "Opaque" : "Transparent", -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_AlphaFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}

	raytracer_app->m_SpecularFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->specularmap.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->specularmap[i].filename);
		wxListItem item;
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->specularmap[i].filename;
		long index = raytracer_app->m_SpecularFile->InsertItem(i, item);
		raytracer_app->m_SpecularFile->SetItem(index, 0, txt, -1);
		raytracer_app->m_SpecularFile->SetItem(index, 1, raytracer.s.objects[selected_object]->specularmap[i].multiplier.toColorStr(), -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_SpecularFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}

	raytracer_app->m_RoughnessFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->roughnessmap.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->roughnessmap[i].filename);
		wxListItem item;
		item.SetText(txt);
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->roughnessmap[i].filename;
		long index = raytracer_app->m_RoughnessFile->InsertItem(i, item);
		raytracer_app->m_RoughnessFile->SetItem(index, 0, txt, -1);
		raytracer_app->m_RoughnessFile->SetItem(index, 1, raytracer.s.objects[selected_object]->roughnessmap[i].multiplier.toRedValueStr(), -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_RoughnessFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}

	raytracer_app->m_TranspFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->transparent_map.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->transparent_map[i].filename);
		wxListItem item;
		item.SetText(txt);
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->transparent_map[i].filename;
		long index = raytracer_app->m_TranspFile->InsertItem(i, item);
		raytracer_app->m_TranspFile->SetItem(index, 0, txt, -1);
		raytracer_app->m_TranspFile->SetItem(index, 1, (raytracer.s.objects[selected_object]->transparent_map[i].multiplier[0]>0.5) ? "Not Refractive" : "Refractive", -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_TranspFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}
	raytracer_app->m_RefrFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->refr_index_map.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->refr_index_map[i].filename);
		wxListItem item;
		item.SetText(txt);
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->refr_index_map[i].filename;
		long index = raytracer_app->m_RefrFile->InsertItem(i, item);
		raytracer_app->m_RefrFile->SetItem(index, 0, txt, -1);
		raytracer_app->m_RefrFile->SetItem(index, 1, raytracer.s.objects[selected_object]->refr_index_map[i].multiplier.toRedValueStr(), -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_RefrFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}

	if (raytracer.s.objects[selected_object]->type == OT_TRIMESH) {
		TriMesh* g = dynamic_cast<TriMesh*>(raytracer.s.objects[selected_object]);
		int id = g->indices[selected_tri].group;
		raytracer_app->m_AlbedoFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		raytracer_app->m_NormalFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		raytracer_app->m_SubsurfaceFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		raytracer_app->m_AlphaFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		raytracer_app->m_SpecularFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		raytracer_app->m_RoughnessFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		raytracer_app->m_RefrFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		raytracer_app->m_TranspFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);


		raytracer_app->m_AlbedoFile->EnsureVisible(id);
		raytracer_app->m_NormalFile->EnsureVisible(id);
		raytracer_app->m_SubsurfaceFile->EnsureVisible(id);
		raytracer_app->m_AlphaFile->EnsureVisible(id);
		raytracer_app->m_SpecularFile->EnsureVisible(id);
		raytracer_app->m_RoughnessFile->EnsureVisible(id);
		raytracer_app->m_RefrFile->EnsureVisible(id);
		raytracer_app->m_TranspFile->EnsureVisible(id);
	}


}

void RenderPanel::render(wxDC& dc)
{
	if (cur_img.size() == 0) return;
	if (raytracer.stopped) return;
	raytracer_app->progressBar->SetValue(raytracer.current_nb_rays / (double)raytracer.nrays*1000.);
	std::ostringstream os;
	os << "Time per ray: " << raytracer.curTimePerFrame / 1000. << " s";
	raytracer_app->infoPerf->SetLabelText(os.str().c_str());

#if 1
	gamma_corrected_image.resize(raytracer.W*raytracer.H * 3);
	extrapolated_image = raytracer.imagedouble;
	/*computed = raytracer.computed;
	computed2 = raytracer.computed;
	double lastR = 0, lastG = 0, lastB = 0;
	int offsetsi[4] = { -1, 0, 0, 1 };
	int offsetsj[4] = { 0, -1, 1, 0 };
	int nbchanged;
	for (int iter = 0; iter < 10; iter++) {   // Voronoi : lent et moche
	nbchanged = 0;
	#pragma omp parallel for
	for (int i = 0; i < raytracer.H; i++) {
	for (int j = 0; j < raytracer.W; j++) {
	if (!computed2[i*raytracer.W + j]) {
	//for (int k = std::max(0, i - 1); k <= std::min(H - 1, i + 1); k++) {
	//	for (int l = std::max(0, j - 1); l <= std::min(W - 1, j + 1); l++) {
	int rand_shuff = rand();
	for (int id=0; id<4; id++) {
	int k = i + offsetsi[(id+ rand_shuff)%4]; if (k < 0 || k>= raytracer.H) continue;
	int l = j + offsetsj[(id + rand_shuff) % 4]; if (l < 0 || l >= raytracer.W) continue;
	if (computed2[k*raytracer.W + l]) {
	extrapolated_image[(i*raytracer.W + j) * 3 + 0] = extrapolated_image[(k*raytracer.W + l) * 3 + 0];
	extrapolated_image[(i*raytracer.W + j) * 3 + 1] = extrapolated_image[(k*raytracer.W + l) * 3 + 1];
	extrapolated_image[(i*raytracer.W + j) * 3 + 2] = extrapolated_image[(k*raytracer.W + l) * 3 + 2];
	computed[i*raytracer.W + j] = true;
	nbchanged++;
	break;
	}
	}
	}
	}
	}
	if (nbchanged == 0) break;
	computed2 = computed;
	}*/

	raytracer.sample_count.resize(raytracer.W*raytracer.H); // just in case.
	for (int i = 0; i < raytracer.H; i++) {
		for (int j = 0; j < raytracer.W; j++) {
			//if (!raytracer.computed[i*raytracer.W + j]) {
			//if (j< raytracer.W/2) {
			if (raytracer.sample_count[i*raytracer.W + j] <= 5) {
				double alpha = raytracer.sample_count[i*raytracer.W + j] / 6.;

				int ix = j / 16;
				double fx = (j / 16.) - ix;

				int iy = i / 16;
				double fy = (i / 16.) - iy;

				for (int k = 0; k < 3; k++) {
					//double lowval = raytracer.imagedouble_lowres[(iy * raytracer.Wlr + ix) * 3 + k];  // nearest neighbor interp.
					double lowval;
					if (ix < raytracer.Wlr - 1 && iy < raytracer.Hlr - 1)
						lowval = (raytracer.imagedouble_lowres[(iy*raytracer.Wlr + ix) * 3 + k] * (1 - fx) + raytracer.imagedouble_lowres[(iy*raytracer.Wlr + ix + 1) * 3 + k] * fx) * (1 - fy)
						+ (raytracer.imagedouble_lowres[((iy + 1)*raytracer.Wlr + ix) * 3 + k] * (1 - fx) + raytracer.imagedouble_lowres[((iy + 1)*raytracer.Wlr + ix + 1) * 3 + k] * fx) * fy;
					else
						lowval = raytracer.imagedouble_lowres[(iy*raytracer.Wlr + ix) * 3 + k];
					extrapolated_image[(i*raytracer.W + j) * 3 + k] = extrapolated_image[(i*raytracer.W + j) * 3 + k] * alpha + (1. - alpha)*lowval;
				}
			}
		}
	}
	//os << "Nbchanged: " << nbchanged;
	//raytracer_app->infoPerf->SetLabelText(os.str().c_str());

#pragma omp parallel for 
	for (int i = 0; i < raytracer.H; i++) {
		for (int j = 0; j < raytracer.W; j++) {
			double normalization = 1. / (raytracer.sample_count[i*raytracer.W + j] + 1.);
			/*gamma_corrected_image[(i*raytracer.W + j) * 3 + 0] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 0] / (raytracer.current_nb_rays + 1)), 1 / 2.2));   // rouge
			gamma_corrected_image[(i*raytracer.W + j) * 3 + 1] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 1] / (raytracer.current_nb_rays + 1)), 1 / 2.2)); // vert
			gamma_corrected_image[(i*raytracer.W + j) * 3 + 2] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 2] / (raytracer.current_nb_rays + 1)), 1 / 2.2)); // bleu*/
			gamma_corrected_image[(i*raytracer.W + j) * 3 + 0] = std::min(255., 255.*fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 0] / 196964.699 * normalization), 1 / raytracer.gamma));   // rouge
			gamma_corrected_image[(i*raytracer.W + j) * 3 + 1] = std::min(255., 255.*fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 1] / 196964.699 * normalization), 1 / raytracer.gamma)); // vert
			gamma_corrected_image[(i*raytracer.W + j) * 3 + 2] = std::min(255., 255.*fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 2] / 196964.699 * normalization), 1 / raytracer.gamma)); // bleu

		}
	}

	static int nbcalls = -1; nbcalls++;
	if (nbcalls == 0 || screenImage.GetWidth() != raytracer.W || screenImage.GetHeight() != raytracer.H) {
		screenImage = wxImage(raytracer.W, raytracer.H, &(gamma_corrected_image[0]), true);
	}

	bmpBuf = wxBitmap(screenImage/*, dc*/);


	double scale_x = (double)displayW / raytracer.W;
	double scale_y = (double)displayH / raytracer.H;
	dc.SetUserScale(scale_x, scale_y);
	dc.DrawBitmap(bmpBuf, 0, 0);
	dc.SetUserScale(1.0, 1.0);
#endif
}


void RaytracerFrame::SaveAs(wxCommandEvent &evt) {
	wxFileDialog
		saveFileDialog(this, _("Save SCN file"), "", "",
			"SCN files (*.scn)|*.scn", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (saveFileDialog.ShowModal() == wxID_CANCEL)
		return;     

	render_panel->raytracer.save_scene(saveFileDialog.GetPath());	
}


void RaytracerFrame::SaveImage(wxCommandEvent &evt) {
	wxFileDialog
		saveFileDialog(this, _("Save Image file"), "", "",
			"Image files (*.*)|*.*", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (saveFileDialog.ShowModal() == wxID_CANCEL)
		return;

	int nbr = render_panel->raytracer.current_nb_rays;
	int W = render_panel->raytracer.W;
	int H = render_panel->raytracer.H;
	std::vector<unsigned char> image(W*H * 3);
#pragma omp parallel for 
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			image[((H - i - 1)*W + j) * 3 + 0] = std::min(255., std::max(0., 255.*std::pow(render_panel->raytracer.imagedouble[((H - i - 1)*W + j) * 3 + 0] / 196964.7 / (nbr + 1.), 1 / render_panel->raytracer.gamma)));   // rouge
			image[((H - i - 1)*W + j) * 3 + 1] = std::min(255., std::max(0., 255.*std::pow(render_panel->raytracer.imagedouble[((H - i - 1)*W + j) * 3 + 1] / 196964.7 / (nbr + 1.), 1 / render_panel->raytracer.gamma))); // vert
			image[((H - i - 1)*W + j) * 3 + 2] = std::min(255., std::max(0., 255.*std::pow(render_panel->raytracer.imagedouble[((H - i - 1)*W + j) * 3 + 2] / 196964.7 / (nbr + 1.), 1 / render_panel->raytracer.gamma))); // bleu
		}
	}
	save_image(saveFileDialog.GetPath().c_str(), &image[0], W, H);
}


void RaytracerFrame::ExportMtl(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	TriMesh* g = dynamic_cast<TriMesh*>(render_panel->raytracer.s.objects[obj_id]);
	if (!g) return;

	wxFileDialog
		saveFileDialog(this, _("Export material"), "", "",
			"MTL files (*.mtl)|*.mtl", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (saveFileDialog.ShowModal() == wxID_CANCEL)
		return;

	g->exportMTL(saveFileDialog.GetPath().c_str());
}

void RaytracerFrame::Open(wxCommandEvent &evt) {
	wxFileDialog
		openFileDialog(this, _("Open SCN file"), "", "",
			"SCN files (*.scn)|*.scn", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
	if (openFileDialog.ShowModal() == wxID_CANCEL)
		return;

	render_panel->stop_render();
	render_panel->raytracer.load_scene(openFileDialog.GetPath());
	render_panel->update_gui();	
	render_panel->start_render();
#ifdef USE_EMBREE
  std::this_thread::sleep_for(std::chrono::milliseconds(100) );
	render_panel->stop_render();
	render_panel->start_render();
#endif
}


void RaytracerFrame::ShowMeshInfo(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	TriMesh* g = dynamic_cast<TriMesh*>(render_panel->raytracer.s.objects[obj_id]);
	if (!g) return;

	int nbE, nbNonM, nbBoundaries;
	int nbC = g->getNbConnected(nbE, nbNonM, nbBoundaries);
	int Euler = (int)g->vertices.size() - nbE + (int)g->indices.size();
	int nbRealEdges, nbOthers, nbTri;
	g->findQuads(nbTri, nbOthers, nbRealEdges);

	std::ostringstream os;
	os << "Object ID: " << obj_id << std::endl;
	os << "Nb Triangles in Triangulated Mesh: " << g->indices.size() << std::endl;
	os << "Nb Polygons in Original Mesh: " << nbTri+ nbOthers<<" (including "<<nbTri<<" triangles and "<<nbOthers<<" other polygons)" << std::endl;	
	os << "Nb Vertices: " << g->vertices.size() << std::endl;
	os << "Nb Edges in Triangulated Mesh: " << nbE << std::endl;
	os << "Nb Edges in Original Mesh: " << nbRealEdges << std::endl;
	os << "Nb Connected Components: " << nbC << std::endl;
	os << "Non Manifold Edges: " << nbNonM << std::endl;
	os << "Boundary edges: " << nbBoundaries << std::endl;
	os << "Euler Number: " << Euler << std::endl;
	os << "Genus: " << (2 - Euler)/2 << std::endl;
	os << " BVH max depth: " << g->bvh_depth << std::endl;
	os << " BVH avg depth: " << g->bvh_avg_depth << std::endl;
	os << " BVH nb nodes: " << g->bvh_nb_nodes << std::endl;
	os << " BVH max triangles per leaf: " << g->max_bvh_triangles << std::endl;
	os << " Nb materials: " << g->textures.size() << std::endl;

	wxMessageBox(os.str().c_str(), "Mesh information");


}

void RaytracerFrame::OnPopupClick(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	size_t item_id = reinterpret_cast<size_t>((static_cast<wxMenu *>(evt.GetEventObject())->GetClientData()));
	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_AlbedoFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);

	switch (evt.GetId()) {
	case ID_ALBEDO_DELETE:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_AlbedoFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			render_panel->raytracer.s.objects[obj_id]->remove_texture(firstSel);
		}
		render_panel->start_render();
		break;	
	case ID_MOVEUP:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_textures(item_id, item_id-1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN:
		if (item_id >= render_panel->raytracer_app->m_AlbedoFile->GetItemCount()-1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_textures(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_REMOVE_TEXTURE:
		if (item_id >= render_panel->raytracer_app->m_AlbedoFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->textures[item_id].clear_texture();
		render_panel->start_render();
		break;
	case ID_ADDWHITE:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->add_col_texture(Vector(col.Red()/255., col.Green() / 255., col.Blue() / 255.));
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_COLOR:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(col.Red() / 255., col.Green() / 255., col.Blue() / 255.), item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();			
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_texture(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_GOLD:
		render_panel->stop_render();
		//http://devernay.free.fr/cours/opengl/materials.html
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.75164, 0.60648, 0.22648), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.628281, 0.555802, 0.366065), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(0.4 * 128, 0.4 * 128, 0.4 * 128), item_id);
		render_panel->start_render();
		break;
	case ID_GOLD_NGAN:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.069, 0.0323, 0.00638), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.0738, 0.0434, 0.0104), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(41.9, 41.9, 41.9), item_id);
		render_panel->start_render();
		break;
	case ID_SILVER:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.50754, 0.50754, 0.50754), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.508273, 0.508273, 0.508273), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(0.4 * 128, 0.4 * 128, 0.4 * 128), item_id);
		render_panel->start_render();
		break;
	case ID_SILVER_NGAN:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.0695, 0.0628, 0.0446), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.0742, 0.0615, 0.0412), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(75, 75, 75), item_id);
		render_panel->start_render();
		break;
	case ID_PEARL:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(1, 0.829, 0.829), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.296648, 0.296648, 0.296648), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(0.088 * 128, 0.088 * 128, 0.088 * 128), item_id);
		render_panel->start_render();
		break;
	case ID_PEARL_NGAN:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.189, 0.146, 0.0861), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.0485, 0.0346, 0.0161), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(27.7, 27.7, 27.7), item_id);
		render_panel->start_render();
		break;
	case ID_WHITE_PLASTIC:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.55, 0.55, 0.55), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.70, 0.70, 0.70), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(0.25 * 128, 0.25 * 128, 0.25 * 128), item_id);
		render_panel->start_render();
		break;
	case ID_WHITE_PLASTIC_NGAN: // actually gray
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.102, 0.0887, 0.0573), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.00699, 0.00566, 0.0036), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(1040, 1040, 1040), item_id);
		render_panel->start_render();
		break;
	case ID_CHROME:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.4, 0.4, 0.4), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.774597, 0.774597, 0.774597), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(0.6 * 128, 0.6 * 128, 0.6 * 128), item_id);
		render_panel->start_render();
		break;
	case ID_CHROME_NGAN:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.00817, 0.0063, 0.00474), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.0213, 0.0151, 0.00766), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(17900, 17900, 17900), item_id);
		render_panel->start_render();
		break;
	case ID_BRONZE:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.714, 0.4284, 0.18144), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.393548, 0.271906, 0.166721), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(0.2 * 128, 0.2 * 128, 0.2 * 128), item_id);
		render_panel->start_render();
		break;
	case ID_BRONZE_NGAN: // alum bronze
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.0864, 0.0597, 0.0302), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.015, 0.00818, 0.00381), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(1290, 1290, 1290), item_id);
		render_panel->start_render();
		break;
	case ID_COPPER:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.7038, 0.27048, 0.0828), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.256777, 0.137622, 0.086014), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(0.1 * 128, 0.1 * 128, 0.1 * 128), item_id);
		render_panel->start_render();
		break;
	case ID_COPPER_NGAN:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(0.0749, 0.0414, 0.027), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(0.0756, 0.0437, 0.0202), item_id);
		render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(33200, 33200, 33200), item_id);
		render_panel->start_render();
		break;
	}
	render_panel->update_gui();
}



void RaytracerFrame::OnPopupClickSpecular(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = reinterpret_cast<size_t>((static_cast<wxMenu *>(evt.GetEventObject())->GetClientData()));
	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_SpecularFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);

	switch (evt.GetId()) {
	case ID_ALBEDO_DELETE_SPECULAR:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_SpecularFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			render_panel->raytracer.s.objects[obj_id]->remove_specular(firstSel);
		}		
		render_panel->start_render();
		break;
	case ID_REMOVE_TEXTURE_SPECULAR:
		if (item_id >= render_panel->raytracer_app->m_SpecularFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->specularmap[item_id].clear_texture();
		render_panel->start_render();
		break;
	case ID_MOVEUP_SPECULAR:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_specular(item_id, item_id - 1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN_SPECULAR:
		if (item_id >= render_panel->raytracer_app->m_SpecularFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_specular(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_ADDWHITE_SPECULAR:
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->add_col_specular(Vector(col.Red() / 255., col.Green() / 255., col.Blue() / 255.));
			render_panel->start_render();			
		}
		break;
	case ID_CHANGE_COLOR_SPECULAR:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(col.Red() / 255., col.Green() / 255., col.Blue() / 255.), item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE_SPECULAR:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_specularmap(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}

void RaytracerFrame::OnPopupClickSubsurface(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = reinterpret_cast<size_t>((static_cast<wxMenu *>(evt.GetEventObject())->GetClientData()));
	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_SubsurfaceFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);

	switch (evt.GetId()) {
	case ID_ALBEDO_DELETE_SUBSURFACE:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_SubsurfaceFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			render_panel->raytracer.s.objects[obj_id]->remove_subsurface(firstSel);
		}
		render_panel->start_render();
		break;
	case ID_REMOVE_TEXTURE_SUBSURFACE:
		if (item_id >= render_panel->raytracer_app->m_SubsurfaceFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->subsurface[item_id].clear_texture();
		render_panel->start_render();
		break;
	case ID_MOVEUP_SUBSURFACE:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_subsurface(item_id, item_id - 1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN_SUBSURFACE:
		if (item_id >= render_panel->raytracer_app->m_SpecularFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_subsurface(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_ADDWHITE_SUBSURFACE:
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->add_col_subsurface(Vector(col.Red() / 255., col.Green() / 255., col.Blue() / 255.));
			render_panel->start_render();
		}
		break;
	case ID_CHANGE_COLOR_SUBSURFACE:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_col_subsurface(Vector(col.Red() / 255., col.Green() / 255., col.Blue() / 255.), item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE_SUBSURFACE:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_subsurface(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}

void RaytracerFrame::OnPopupClickRoughness(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = reinterpret_cast<size_t>(static_cast<wxMenu *>(evt.GetEventObject())->GetClientData());
	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_RoughnessFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);

	switch (evt.GetId()) {
	case ID_DELETE_ROUGHNESS:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_RoughnessFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			render_panel->raytracer.s.objects[obj_id]->remove_roughness(firstSel);
		}		
		render_panel->start_render();
		break;
	case ID_REMOVE_TEXTURE_ROUGHNESS:
		if (item_id >= render_panel->raytracer_app->m_RoughnessFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->roughnessmap[item_id].clear_texture();
		render_panel->start_render();
		break;
	case ID_MOVEUP_ROUGHNESS:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_roughness(item_id, item_id - 1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN_ROUGHNESS:
		if (item_id >= render_panel->raytracer_app->m_RoughnessFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_roughness(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_ADDWHITE_ROUGHNESS:
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->add_col_roughness(Vector(col.Red(), col.Green(), col.Blue()));
			render_panel->start_render();			
		}
		break;
	case ID_CHANGE_COLOR_ROUGHNESS:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(col.Red(), col.Green(), col.Blue()), item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE_ROUGHNESS:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_roughnessmap(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}



void RaytracerFrame::OnPopupClickRefrIndex(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = reinterpret_cast<size_t>(static_cast<wxMenu *>(evt.GetEventObject())->GetClientData());
	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_RefrFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
	std::string ret;

	switch (evt.GetId()) {
	case ID_DELETE_REFR:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_RefrFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			render_panel->raytracer.s.objects[obj_id]->remove_refr(firstSel);
		}
		render_panel->start_render();
		break;
	case ID_REMOVE_TEXTURE_REFR:
		if (item_id >= render_panel->raytracer_app->m_RefrFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->refr_index_map[item_id].clear_texture();
		render_panel->start_render();
		break;
	case ID_MOVEUP_REFR:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_refr(item_id, item_id - 1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN_REFR:
		if (item_id >= render_panel->raytracer_app->m_RefrFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_refr(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_ADDWHITE_REFR:
		ret = wxGetTextFromUser("Enter refraction index", "Refraction index", "1.3").ToStdString();
		if (ret.size() != 0) {
			double val;
			try {
				val = std::stod(ret.c_str());
			} catch (...) {
				return;
			}
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->add_col_refr(val);
			render_panel->start_render();
		}
		break;
	case ID_CHANGE_COLOR_REFR:
	{
		//render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK
		ret = wxGetTextFromUser("Enter refraction index", "Refraction index", "1.3");
		if (ret.size()!=0 ) {			
			double val;
			try {
				val = std::stod(ret.c_str());
			} catch (...) {
				return;
			}
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_col_refr(val, item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE_REFR:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_refr_map(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}


void RaytracerFrame::OnPopupClickTransparent(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = reinterpret_cast<size_t>(static_cast<wxMenu *>(evt.GetEventObject())->GetClientData());
	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_TranspFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
	wxArrayString choices; choices.push_back("Yes"); choices.push_back("No");

	switch (evt.GetId()) {
	case ID_DELETE_TRANSP:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_TranspFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			render_panel->raytracer.s.objects[obj_id]->remove_transp(firstSel);
		}
		render_panel->start_render();
		break;
	case ID_REMOVE_TEXTURE_TRANSP:
		if (item_id >= render_panel->raytracer_app->m_TranspFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->transparent_map[item_id].clear_texture();
		render_panel->start_render();
		break;
	case ID_MOVEUP_TRANSP:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_transp(item_id, item_id - 1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN_TRANSP:
		if (item_id >= render_panel->raytracer_app->m_TranspFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_transp(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_ADDWHITE_TRANSP:
	{		
		int ret = wxGetSingleChoiceIndex("Is this material at least partly transparent ?", "Transparent ?", choices, 1);
		if (ret!=-1) {
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->add_col_transp(ret==0?0:1);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_COLOR_TRANSP:
	{
		int ret = wxGetSingleChoiceIndex("Is this material at least partly transparent ?", "Transparent ?", choices, 1);
		if (ret != -1) {
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_col_transp(ret == 0 ? 0 : 1, item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE_TRANSP:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_transp_map(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}

void RaytracerFrame::OnPopupClickNormal(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = reinterpret_cast<size_t>(static_cast<wxMenu *>(evt.GetEventObject())->GetClientData());

	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_NormalFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
	switch (evt.GetId()) {
	case ID_ALBEDO_DELETE_NORMAL:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_NormalFile->GetNextItem(itemIndex, 	wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			//wxLogDebug(listControl->GetItemText(itemIndex));
			render_panel->raytracer.s.objects[obj_id]->remove_normal(firstSel);
		}
		//render_panel->raytracer.s.objects[obj_id]->remove_normal(item_id);
		render_panel->start_render();
		break;
	case ID_MOVEUP_NORMAL:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_normal(item_id, item_id - 1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN_NORMAL:
		if (item_id >= render_panel->raytracer_app->m_NormalFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_normal(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_ADDNULL_NORMAL:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->add_null_normalmap();
		render_panel->start_render();
		break;
	case ID_CHANGE_COLOR_NORMAL:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_null_normalmap(item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE_NORMAL:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_normalmap(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}
void RaytracerFrame::OnPopupClickAlpha(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = reinterpret_cast<size_t>(static_cast<wxMenu *>(evt.GetEventObject())->GetClientData());
	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_AlphaFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
	switch (evt.GetId()) {
	case ID_ALBEDO_DELETE_ALPHA:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_AlphaFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			render_panel->raytracer.s.objects[obj_id]->remove_alpha(firstSel);
		}
		render_panel->start_render();
		break;
	case ID_MOVEUP_ALPHA:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_alpha(item_id, item_id - 1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN_ALPHA:
		if (item_id >= render_panel->raytracer_app->m_AlphaFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_alpha(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_ADDWHITE_ALPHA:
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->add_col_alpha(col.Red() / 255.);
			render_panel->start_render();		
		}
		break;
	case ID_CHANGE_COLOR_ALPHA:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_col_alpha(col.Red() / 255., item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE_ALPHA:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_alphamap(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}

void RaytracerFrame::DeselectAll(wxListCtrl* list) {

	long item = -1;
	for (;; ) {
		item = list->GetNextItem(item, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
		if (item == -1)
			break;

		list->SetItemState(item, 0, wxLIST_STATE_SELECTED);
	}
}

void RaytracerFrame::OnListSelected(wxListEvent &evt) {
	//int sel = evt.GetItem();
	std::vector<int> selected_items;
	long item = -1;
	for (;; ) {
		item = ((wxListCtrl*)evt.GetEventObject())->GetNextItem(item, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
		if (item == -1)
			break;

		selected_items.push_back(item);
	}
	
	if (programHandling) { return; }  // this is because SetItemState invokes an OnSelectedItem event, which in turn calls this function indefinitely
	programHandling = true;
	DeselectAll(render_panel->raytracer_app->m_AlbedoFile);
	DeselectAll(render_panel->raytracer_app->m_NormalFile);
	DeselectAll(render_panel->raytracer_app->m_AlphaFile);
	DeselectAll(render_panel->raytracer_app->m_SpecularFile);
	DeselectAll(render_panel->raytracer_app->m_RoughnessFile);
	DeselectAll(render_panel->raytracer_app->m_RefrFile);
	DeselectAll(render_panel->raytracer_app->m_TranspFile);

	for (int i = 0; i < selected_items.size(); i++) {
		render_panel->raytracer_app->m_AlbedoFile->SetItemState(selected_items[i], wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		render_panel->raytracer_app->m_NormalFile->SetItemState(selected_items[i], wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		render_panel->raytracer_app->m_AlphaFile->SetItemState(selected_items[i], wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		render_panel->raytracer_app->m_SpecularFile->SetItemState(selected_items[i], wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		render_panel->raytracer_app->m_RoughnessFile->SetItemState(selected_items[i], wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		render_panel->raytracer_app->m_RefrFile->SetItemState(selected_items[i], wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		render_panel->raytracer_app->m_TranspFile->SetItemState(selected_items[i], wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
	}

	int id = evt.GetItem();
	render_panel->raytracer_app->m_AlbedoFile->EnsureVisible(id);
	render_panel->raytracer_app->m_NormalFile->EnsureVisible(id);
	render_panel->raytracer_app->m_AlphaFile->EnsureVisible(id);
	render_panel->raytracer_app->m_SpecularFile->EnsureVisible(id);
	render_panel->raytracer_app->m_RoughnessFile->EnsureVisible(id);
	render_panel->raytracer_app->m_RefrFile->EnsureVisible(id);
	render_panel->raytracer_app->m_TranspFile->EnsureVisible(id);
	evt.Skip();
	programHandling = false;
}

void RaytracerFrame::createAlbedoMenu() {
	
	albedomnu.Append(ID_MOVEUP, "Move Up");
	albedomnu.Append(ID_MOVEDOWN, "Move Down");
	albedomnu.Append(ID_ALBEDO_DELETE, "Delete Row");
	albedomnu.Append(ID_REMOVE_TEXTURE, "Remove Texture");
	albedomnu.Append(ID_ADDWHITE, "Add Color");
	albedomnu.Append(ID_CHANGE_COLOR, "Change Color");
	albedomnu.Append(ID_CHANGE_TEXTURE, "Change Texture");

	subsubmnu = new wxMenu();
	albedosubmnu = new wxMenu();
	subsubmnu2 = new wxMenu();
	subsubmnu3 = new wxMenu();
	subsubmnu4 = new wxMenu();
	subsubmnu5 = new wxMenu();
	subsubmnu6 = new wxMenu();
	subsubmnu7 = new wxMenu();

	subsubmnu->Append(ID_GOLD_NGAN, "Ngan");
	subsubmnu->Append(ID_GOLD, "OpenGL");
	albedosubmnu->AppendSubMenu(subsubmnu, "Gold");

	subsubmnu2->Append(ID_SILVER_NGAN, "Ngan");
	subsubmnu2->Append(ID_SILVER, "OpenGL");
	albedosubmnu->AppendSubMenu(subsubmnu2, "Silver");

	subsubmnu3->Append(ID_CHROME_NGAN, "Ngan");
	subsubmnu3->Append(ID_CHROME, "OpenGL");
	albedosubmnu->AppendSubMenu(subsubmnu3, "Chrome");

	subsubmnu4->Append(ID_BRONZE_NGAN, "Ngan");
	subsubmnu4->Append(ID_BRONZE, "OpenGL");
	albedosubmnu->AppendSubMenu(subsubmnu4, "Bronze");

	subsubmnu5->Append(ID_COPPER_NGAN, "Ngan");
	subsubmnu5->Append(ID_COPPER, "OpenGL");
	albedosubmnu->AppendSubMenu(subsubmnu5, "Copper");

	subsubmnu6->Append(ID_WHITE_PLASTIC_NGAN, "Ngan");
	subsubmnu6->Append(ID_WHITE_PLASTIC, "OpenGL");
	albedosubmnu->AppendSubMenu(subsubmnu6, "White Plastic");

	subsubmnu7->Append(ID_PEARL_NGAN, "Ngan");
	subsubmnu7->Append(ID_PEARL, "OpenGL");
	albedosubmnu->AppendSubMenu(subsubmnu7, "Pearl");

	albedomnu.AppendSubMenu(albedosubmnu, "Change Whole Material Color to...");
	albedomnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClick), NULL, this);
}

void RaytracerFrame::OnListRightClick(wxListEvent &evt) {
	size_t sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	albedomnu.SetClientData(data);
	PopupMenu(&albedomnu);
}

void RaytracerFrame::OnListRightClickSpecular(wxListEvent &evt) {
	size_t sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP_SPECULAR, "Move Up");
	mnu.Append(ID_MOVEDOWN_SPECULAR, "Move Down");
	mnu.Append(ID_ALBEDO_DELETE_SPECULAR, "Delete Row");
	mnu.Append(ID_REMOVE_TEXTURE_SPECULAR, "Remove Texture");
	mnu.Append(ID_ADDWHITE_SPECULAR, "Add Color");
	mnu.Append(ID_CHANGE_COLOR_SPECULAR, "Change Color");
	mnu.Append(ID_CHANGE_TEXTURE_SPECULAR, "Change Texture");

	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClickSpecular), NULL, this);
	PopupMenu(&mnu);
}

void RaytracerFrame::OnListRightClickSubsurface(wxListEvent &evt) {
	size_t sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP_SUBSURFACE, "Move Up");
	mnu.Append(ID_MOVEDOWN_SUBSURFACE, "Move Down");
	mnu.Append(ID_ALBEDO_DELETE_SUBSURFACE, "Delete Row");
	mnu.Append(ID_REMOVE_TEXTURE_SUBSURFACE, "Remove Texture");
	mnu.Append(ID_ADDWHITE_SUBSURFACE, "Add Color");
	mnu.Append(ID_CHANGE_COLOR_SUBSURFACE, "Change Color");
	mnu.Append(ID_CHANGE_TEXTURE_SUBSURFACE, "Change Texture");

	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClickSubsurface), NULL, this);
	PopupMenu(&mnu);
}
void RaytracerFrame::OnListRightClickNormal(wxListEvent &evt) {
	size_t sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP_NORMAL, "Move Up");
	mnu.Append(ID_MOVEDOWN_NORMAL, "Move Down");
	mnu.Append(ID_ALBEDO_DELETE_NORMAL, "Delete");
	mnu.Append(ID_ADDNULL_NORMAL, "Add Null");
	mnu.Append(ID_CHANGE_COLOR_NORMAL, "Change to Null");
	mnu.Append(ID_CHANGE_TEXTURE_NORMAL, "Change Texture");

	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClickNormal), NULL, this);
	PopupMenu(&mnu);
}
void RaytracerFrame::OnListRightClickAlpha(wxListEvent &evt) {
	size_t sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP_ALPHA, "Move Up");
	mnu.Append(ID_MOVEDOWN_ALPHA, "Move Down");
	mnu.Append(ID_ALBEDO_DELETE_ALPHA, "Delete");
	mnu.Append(ID_ADDWHITE_ALPHA, "Add Color");
	mnu.Append(ID_CHANGE_COLOR_ALPHA, "Change Color");
	mnu.Append(ID_CHANGE_TEXTURE_ALPHA, "Change Texture");

	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClickAlpha), NULL, this);
	PopupMenu(&mnu);
}
void RaytracerFrame::OnListRightClickRoughness(wxListEvent &evt) {
	size_t sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP_ROUGHNESS, "Move Up");
	mnu.Append(ID_MOVEDOWN_ROUGHNESS, "Move Down");
	mnu.Append(ID_DELETE_ROUGHNESS, "Delete Row");
	mnu.Append(ID_REMOVE_TEXTURE_ROUGHNESS, "Remove Texture");
	mnu.Append(ID_ADDWHITE_ROUGHNESS, "Add Color");
	mnu.Append(ID_CHANGE_COLOR_ROUGHNESS, "Change Color");
	mnu.Append(ID_CHANGE_TEXTURE_ROUGHNESS, "Change Texture");
	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClickRoughness), NULL, this);
	PopupMenu(&mnu);
}

void RaytracerFrame::OnListRightClickRefrIndex(wxListEvent &evt) {
	size_t sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP_REFR, "Move Up");
	mnu.Append(ID_MOVEDOWN_REFR, "Move Down");
	mnu.Append(ID_DELETE_REFR, "Delete Row");
	mnu.Append(ID_REMOVE_TEXTURE_REFR, "Remove Texture");
	mnu.Append(ID_ADDWHITE_REFR, "Add Index");
	mnu.Append(ID_CHANGE_COLOR_REFR, "Change Index");
	mnu.Append(ID_CHANGE_TEXTURE_REFR, "Change Texture");
	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClickRefrIndex), NULL, this);
	PopupMenu(&mnu);
}
void RaytracerFrame::OnListRightClickTransparent(wxListEvent &evt) {
	size_t sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP_TRANSP, "Move Up");
	mnu.Append(ID_MOVEDOWN_TRANSP, "Move Down");
	mnu.Append(ID_DELETE_TRANSP, "Delete Row");
	mnu.Append(ID_REMOVE_TEXTURE_TRANSP, "Remove Texture");
	mnu.Append(ID_ADDWHITE_TRANSP, "Add Color");
	mnu.Append(ID_CHANGE_COLOR_TRANSP, "Change Color");
	mnu.Append(ID_CHANGE_TEXTURE_TRANSP, "Change Texture");
	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClickTransparent), NULL, this);
	PopupMenu(&mnu);
}
bool DnDFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_pOwner != NULL) {
		m_pOwner->render_panel->stop_render();
		bool has_loaded_obj = false;
		for (size_t n = 0; n < nFiles; n++) {			

			if (filenames[n].Lower().find(".seg") != std::string::npos) {
				if (m_pOwner->render_panel->selected_object < 0 || m_pOwner->render_panel->selected_object>= m_pOwner->render_panel->raytracer.s.objects.size()) {
					wxMessageBox("No triangle mesh was selected to apply the segments to", "Error", wxOK) ;
					return true;
				}
				TriMesh* g = dynamic_cast<TriMesh*>(m_pOwner->render_panel->raytracer.s.objects[m_pOwner->render_panel->selected_object]);
				if (!g) {
					wxMessageBox("Selected object is not a triangle mesh", "Error", wxOK);
					return true;
				}

				std::vector<int> unpermuted_index(g->indices.size());
				for (int i = 0; i < g->indices.size(); i++) {
					unpermuted_index[g->permuted_triangle_index[i]] = i;
				}
				FILE* f = fopen(filenames[n].c_str(), "r+");
				int u;
				int faceid = 0;
				g->facecolors.resize(g->indices.size());
				while (fscanf(f, "%u", &u)!= EOF) {
					Vector col(((u*u*(u + 2) * 123 + 51) % 1000) / 1000., ((u*(u+7)*456 + 266) % 1000) / 1000., ((u*u*u*5+u*33 + 687) % 1000) / 1000.);
					if (faceid<g->indices.size())
						g->facecolors[unpermuted_index[faceid]] = col;
					faceid++;
				}
				fclose(f);
				continue;
			}

			if (filenames[n].Lower().find(".lab") != std::string::npos) {
				if (m_pOwner->render_panel->selected_object < 0 || m_pOwner->render_panel->selected_object >= m_pOwner->render_panel->raytracer.s.objects.size()) {
					wxMessageBox("No triangle mesh was selected to apply the segments to", "Error", wxOK);
					return true;
				}
				TriMesh* g = dynamic_cast<TriMesh*>(m_pOwner->render_panel->raytracer.s.objects[m_pOwner->render_panel->selected_object]);
				if (!g) {
					wxMessageBox("Selected object is not a triangle mesh", "Error", wxOK);
					return true;
				}

				std::vector<int> unpermuted_index(g->indices.size());
				for (int i = 0; i < g->indices.size(); i++) {
					unpermuted_index[g->permuted_triangle_index[i]] = i;
				}
				FILE* f = fopen(filenames[n].c_str(), "r+");
				int u;
				int faceid = 0;
				g->facecolors.resize(g->indices.size());
				char lineName[512];
				char lineFaces[512*1000];
				char* line;
				int segId = 0;
				while (fscanf(f, "%[^\n]\n%[^\n]\n", lineName, lineFaces) != EOF) {
					int faceid, offset;
					line = lineFaces;
					while (sscanf(line, "%u%n", &faceid,&offset) == 1) {						
						faceid--;
						line = line + offset;
						Vector col(((segId*segId*(segId + 2) * 123 + 51) % 1000) / 1000., ((segId*(segId + 7) * 456 + 266) % 1000) / 1000., ((segId*segId*segId * 5 + segId * 33 + 687) % 1000) / 1000.);
						if (faceid < g->indices.size())
							g->facecolors[unpermuted_index[faceid]] = col;
					}
					segId++;
				}
				fclose(f);
				continue;
			}

			if (filenames[n].Lower().find(".xyz") != std::string::npos) {

				std::string values = wxGetTextFromUser("Space separated, -1: ignore, 0:x 1:y, 2:z, 3:nx,4:ny,5:nz, 6:colr, 7:colg, 8:colb",
					"Enter XYZ file format", "-1 -1 0 1 2 6 7 8").ToStdString();
				bool center = (wxMessageBox("Should the model be normalized / centered ?", "Normalization ?", wxYES_NO | wxCENTRE) == wxYES);
				const char* line = &values.c_str()[0];
				int cols[100];
				int nbcols = 0;
				int offset;
				while (sscanf(line, "%d%n", &cols[nbcols], &offset) == 1) {
					nbcols++;
					line += offset;
				}

				PointSet* g = new PointSet(filenames[n], nbcols, cols, false, false, center);
				g->scale = 30;
				Vector c = (g->bvh.bbox.bounds[0] + g->bvh.bbox.bounds[1])*0.5; // c is at 0
				double s = std::max(g->bvh.bbox.bounds[1][0] - g->bvh.bbox.bounds[0][0], std::max(g->bvh.bbox.bounds[1][1] - g->bvh.bbox.bounds[0][1], g->bvh.bbox.bounds[1][2] - g->bvh.bbox.bounds[0][2]));
				g->max_translation = Vector(0, m_pOwner->render_panel->raytracer.s.objects[2]->get_translation(m_pOwner->render_panel->raytracer.s.current_time, m_pOwner->render_panel->raytracer.is_recording)[1] - (g->bvh.bbox.bounds[0][1])*g->scale, 0);
				m_pOwner->render_panel->raytracer.s.addObject(g);
				continue;
			} 

			if ((filenames[n].Lower().find(".obj") != std::string::npos)  || (filenames[n].Lower().find(".wrl") != std::string::npos) || (filenames[n].Lower().find(".off") != std::string::npos)) {
				bool center = (wxMessageBox("Should the model be normalized / centered ?", "Normalization ?", wxYES_NO | wxCENTRE) == wxYES);
				TriMesh* g = new TriMesh(&m_pOwner->render_panel->raytracer.s, filenames[n], 1, Vector(0, 0, 0), false, NULL, false, center);
				g->scale = 30;
				g->display_edges = false;
				has_loaded_obj = true;
				g->max_translation = Vector(0, m_pOwner->render_panel->raytracer.s.objects[2]->get_translation(m_pOwner->render_panel->raytracer.s.current_time, m_pOwner->render_panel->raytracer.is_recording)[1] - (g->bbox.bounds[0][1])*g->scale, 0);
				m_pOwner->render_panel->raytracer.s.addObject(g);
				continue;
			}

			if (filenames[n].Lower().find(".yarn") != std::string::npos) {
				Yarns* y = new Yarns(filenames[n]);
				m_pOwner->render_panel->raytracer.s.addObject(y);
			}

			if (filenames[n].Lower().find(".titopo") != std::string::npos) {
				if (m_pOwner->render_panel->selected_object < 0 || m_pOwner->render_panel->selected_object >= m_pOwner->render_panel->raytracer.s.objects.size()) {
					wxMessageBox("No object was selected to apply the BRDF to", "Error", wxOK);
					return true;
				}
				if (filenames[n].Lower().find(".titopoh") != std::string::npos) {
					m_pOwner->render_panel->raytracer.s.objects[m_pOwner->render_panel->selected_object]->brdf = new TitopoBRDF(filenames[n].ToStdString(), 45, 45, 180);
				} else {
					m_pOwner->render_panel->raytracer.s.objects[m_pOwner->render_panel->selected_object]->brdf = new TitopoBRDF(filenames[n].ToStdString(), 90, 90, 360);
				}
				continue;
			}
			if (filenames[n].Lower().find(".binary") != std::string::npos) {
				if (m_pOwner->render_panel->selected_object < 0 || m_pOwner->render_panel->selected_object >= m_pOwner->render_panel->raytracer.s.objects.size()) {
					wxMessageBox("No object was selected to apply the BRDF to", "Error", wxOK);
					return true;
				}

				m_pOwner->render_panel->raytracer.s.objects[m_pOwner->render_panel->selected_object]->brdf = new IsoMERLBRDF(filenames[n].ToStdString());
				continue;
			}
		}		
		m_pOwner->render_panel->start_render();

	}
	return true;
}

bool DnDSpecularFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_specularmap(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}

bool DnDSubsurfaceFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames) {
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_subsurface(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}


bool DnDAlbedoFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_texture(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}

bool DnDNormalFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_normalmap(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}
bool DnDRoughnessFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_roughnessmap(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}
bool DnDRefrFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_refr_map(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}
bool DnDAlphaFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_alphamap(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}

bool DnDTranspFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_transp_map(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}
bool DnDEnvmapFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (nFiles < 1) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		((Sphere*)m_rtFrame->render_panel->raytracer.s.objects[1])->load_envmap(filenames[0]);
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}

bool DnDBackgroundFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (nFiles < 1) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		m_rtFrame->render_panel->raytracer.s.load_background(filenames[0], m_rtFrame->render_panel->raytracer.gamma);
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}
void RaytracerFrame::OnClose(wxCloseEvent& event)
{
	render_panel->stop_render();
	Destroy();  // you may also do:  event.Skip();
				// since the default event handler does call Destroy(), too
}
