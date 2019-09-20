
#include "wx/wx.h"
#include "wx/dnd.h"
#include "wx/dirctrl.h"
#include "wx/notebook.h"
#include "wx/listctrl.h"
#include <wx/gauge.h>
#include <wx/clrpicker.h>
#include <wx/colordlg.h>
#include <wx/spinctrl.h>
#include <wx/button.h>
#include <wx/tglbtn.h>
#include <ostream>
#include <sstream>
#include "chrono.h"
//#include <Windows.h>
#include <thread>

#include "Raytracer.h"

#define EDGES_CHECKBOX 1001
#define INTERP_CHECKBOX 1002
#define FOV_SLIDER 1003
#define APERTURE_SLIDER 1004
#define KS_SLIDER 1005
#define ALBEDO_FILES 1006
#define SPECULAR_FILES 10061
#define NORMAL_FILES 10062
#define ALPHA_FILES 10063
#define ROUGHNESS_FILES 10064
#define TRANSP_FILES 10065
#define REFR_FILES 10066
#define FILTER_SLIDER 1007
#define ALBEDO_COLORPICKER 1008
#define FOGDENSITY_SLIDER 1009
#define UNIFORMFOG_RADIO 1010
#define EXPFOG_RADIO 1010
#define ENVMAPINTENSITY_SLIDER 1011
#define LIGHTINTENSITY_SLIDER 1012
#define BOUNCES_SPIN 1013
#define TRANSPARENT_CHECKBOX 1014
#define FLIPNORMALS_CHECKBOX 1015
#define GHOST_CHECKBOX 10151
#define REFRACTION_INDEX 1016
#define DELETE_OBJECT 1017
#define RENDER_WIDTH 1018
#define RENDER_HEIGHT 1019
#define FOCUS_SLIDER 1020
#define NBRAYS 1021
#define LAUNCH_RENDER 1022
#define IS_LENTICULAR_CHECKBOX 1023
#define MAXANGLE_SLIDER 1024
#define NBVIEWS 1025
#define NBPIXPERSLICE 1026
#define NBFRAMES 1027
#define ADD_KEYFRAME 1028
#define TIME_SLIDER 1029
#define DURATION 1030
#define RECORD_KEYFRAME 1031
#define RENDER_VIDEO 1032
#define COLOR_ANISOTROPY 1033
#define RANDOM_COLOR 1034

#define ID_ALBEDO_DELETE 10
#define ID_MOVEUP 11
#define ID_MOVEDOWN 12
#define ID_ADDWHITE 13
#define ID_ALBEDO_DELETE_ALPHA 14
#define ID_MOVEUP_ALPHA 15
#define ID_MOVEDOWN_ALPHA 16
#define ID_ADDWHITE_ALPHA 17
#define ID_ALBEDO_DELETE_NORMAL 18
#define ID_MOVEUP_NORMAL 19
#define ID_MOVEDOWN_NORMAL 20
#define ID_ADDNULL_NORMAL 21
#define ID_ALBEDO_DELETE_SPECULAR 22
#define ID_MOVEUP_SPECULAR 23
#define ID_MOVEDOWN_SPECULAR 24
#define ID_ADDWHITE_SPECULAR 25
#define ID_DELETE_ROUGHNESS 26
#define ID_MOVEUP_ROUGHNESS 27
#define ID_MOVEDOWN_ROUGHNESS 28
#define ID_ADDWHITE_ROUGHNESS 29
#define ID_DELETE_TRANSP 261
#define ID_MOVEUP_TRANSP 271
#define ID_MOVEDOWN_TRANSP 281
#define ID_ADDWHITE_TRANSP 291
#define ID_REMOVE_TEXTURE_TRANSP 431
#define ID_CHANGE_COLOR_TRANSP 401
#define ID_CHANGE_TEXTURE_TRANSP 411
#define ID_DELETE_REFR 262
#define ID_MOVEUP_REFR 272
#define ID_MOVEDOWN_REFR 282
#define ID_ADDWHITE_REFR 292
#define ID_REMOVE_TEXTURE_REFR 432
#define ID_CHANGE_COLOR_REFR 402
#define ID_CHANGE_TEXTURE_REFR 412

#define ID_SAVE 30
#define ID_OPEN 31
#define ID_CHANGE_COLOR 32
#define ID_CHANGE_TEXTURE 33
#define ID_CHANGE_COLOR_SPECULAR 34
#define ID_CHANGE_TEXTURE_SPECULAR 35
#define ID_CHANGE_COLOR_NORMAL 36
#define ID_CHANGE_TEXTURE_NORMAL 37
#define ID_CHANGE_COLOR_ALPHA 38
#define ID_CHANGE_TEXTURE_ALPHA 39
#define ID_CHANGE_COLOR_ROUGHNESS 40
#define ID_CHANGE_TEXTURE_ROUGHNESS 41
#define ID_REMOVE_TEXTURE 42
#define ID_REMOVE_TEXTURE_ROUGHNESS 43
#define ID_REMOVE_TEXTURE_SPECULAR 44
#define ID_SAVEIMAGE 45
#define ID_EXPORT_MTL 46
#define ID_MESHINFO 47

#define ID_GOLD 49
#define ID_GOLD_NGAN 491
#define ID_SILVER 50
#define ID_SILVER_NGAN 501
#define ID_PEARL 51
#define ID_PEARL_NGAN 511
#define ID_WHITE_PLASTIC 52
#define ID_WHITE_PLASTIC_NGAN 521
#define ID_CHROME 53
#define ID_CHROME_NGAN 531
#define ID_BRONZE 54
#define ID_BRONZE_NGAN 541
#define ID_COPPER 55
#define ID_COPPER_NGAN 551


#if !defined(wxOVERRIDE)
		#define wxOVERRIDE override;
#endif

class RaytracerApp;
class RenderPanel;

class RaytracerFrame : public wxFrame
{
public:
	RaytracerFrame();
	virtual ~RaytracerFrame() {};


	void OnPaint(wxPaintEvent& event) {};
	void OnSize(wxSizeEvent& event) {
		Refresh();
		event.Skip();
	}
	void OnClose(wxCloseEvent& event);

	void OnLeftDown(wxMouseEvent& event) {};
	void OnRightDown(wxMouseEvent& event) {};

	void OnUpdateUIMoveByDefault(wxUpdateUIEvent& event) {};

	void OnPopupClick(wxCommandEvent &evt);
	void OnPopupClickSpecular(wxCommandEvent &evt);
	void OnPopupClickNormal(wxCommandEvent &evt);
	void OnPopupClickAlpha(wxCommandEvent &evt);
	void OnPopupClickRoughness(wxCommandEvent &evt);
	void OnPopupClickTransparent(wxCommandEvent &evt);
	void OnPopupClickRefrIndex(wxCommandEvent &evt);
	void OnListRightClick(wxListEvent &evt);
	void OnListRightClickSpecular(wxListEvent &evt);
	void OnListRightClickNormal(wxListEvent &evt);
	void OnListRightClickAlpha(wxListEvent &evt);
	void OnListRightClickRoughness(wxListEvent &evt);
	void OnListRightClickTransparent(wxListEvent &evt);
	void OnListRightClickRefrIndex(wxListEvent &evt);
	void OnListSelected(wxListEvent &evt);


	void SaveAs(wxCommandEvent &evt);
	void SaveImage(wxCommandEvent &evt);
	void ExportMtl(wxCommandEvent &evt);
	void Open(wxCommandEvent &evt);
	void ShowMeshInfo(wxCommandEvent &evt);
	void DeselectAll(wxListCtrl* list);

	RenderPanel *render_panel;
	bool programHandling;
private:
	void createAlbedoMenu();
	wxMenu albedomnu, *albedosubmnu;
	wxMenu *subsubmnu, *subsubmnu2, *subsubmnu3, *subsubmnu4, *subsubmnu5, *subsubmnu6, *subsubmnu7;

	wxDECLARE_EVENT_TABLE();
};

class DnDFile : public wxFileDropTarget
{
public:
    DnDFile(RaytracerFrame *pOwner = NULL) { m_pOwner = pOwner; }

    virtual bool OnDropFiles(wxCoord x, wxCoord y,
                             const wxArrayString& filenames) ;

private:
	RaytracerFrame *m_pOwner;
};



class DnDAlbedoFile : public wxFileDropTarget
{
public:
	DnDAlbedoFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) ;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

class DnDSpecularFile : public wxFileDropTarget
{
public:
	DnDSpecularFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) ;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

class DnDRoughnessFile : public wxFileDropTarget
{
public:
	DnDRoughnessFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) ;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};
class DnDRefrFile : public wxFileDropTarget
{
public:
	DnDRefrFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) wxOVERRIDE;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};
class DnDEnvmapFile : public wxFileDropTarget
{
public:
	DnDEnvmapFile(wxTextCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) ;

private:
	wxTextCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

class DnDBackgroundFile : public wxFileDropTarget
{
public:
	DnDBackgroundFile(wxTextCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames);

private:
	wxTextCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};


class DnDNormalFile : public wxFileDropTarget
{
public:
	DnDNormalFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) ;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

class DnDAlphaFile : public wxFileDropTarget
{
public:
	DnDAlphaFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) ;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

class DnDTranspFile : public wxFileDropTarget
{
public:
	DnDTranspFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) wxOVERRIDE;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

wxBEGIN_EVENT_TABLE(RaytracerFrame, wxFrame)
EVT_LEFT_DOWN(RaytracerFrame::OnLeftDown)
EVT_RIGHT_DOWN(RaytracerFrame::OnRightDown)
EVT_PAINT(RaytracerFrame::OnPaint)
EVT_CLOSE(RaytracerFrame::OnClose)
EVT_SIZE(RaytracerFrame::OnSize)
wxEND_EVENT_TABLE()



static void doRender(void* Param);



class RenderPanel : public wxPanel
{

public:

	RenderPanel(RaytracerFrame* parent, RaytracerApp* app):wxPanel(parent) {
		Connect(wxEVT_PAINT, wxPaintEventHandler(RenderPanel::paintEvent));

		raytracer_app = app;
		should_terminate = false;
		raytracer.loadScene();
		W = raytracer.W;
		H = raytracer.H;
		cur_img.resize(W*H * 3);
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++) {
				cur_img[(i*W + j) * 3 + 0] = 255;
				cur_img[(i*W + j) * 3 + 1] = 0;
				cur_img[(i*W + j) * 3 + 2] = 0;
			}
		}

		left_mouse_down = false;
		parent->render_panel = this;
		wasDragging = false;

		//compute_thread = CreateThread(NULL, 0, doRender, (void*) this, 0, NULL);
		//SetThreadPriority(compute_thread, THREAD_PRIORITY_HIGHEST);
		start_render();
		//SetThreadPriority(compute_thread, THREAD_PRIORITY_LOWEST);
	}


	/*void partition(int* indices, int i0, int i1) {

		int mini = 100, maxi = -100;
		for (int i = i0; i < i1; i++) {
			mini = std::min(mini, indices[i]);
			maxi = std::max(maxi, indices[i]);
		}
		double split_val = (mini + maxi) / 2;

		int pivot = i0 - 1;

		for (int i = i0; i < i1; i++) {
			double center_split_dim = indices[i];

			if (center_split_dim <= split_val) {
				pivot++;
				std::swap(indices[i], indices[pivot]);
			}
		}

		if (pivot < i0 || pivot >= i1 - 1 || i1 <= i0 + 1) {
			return;
		}

		partition(indices, i0, pivot + 1);
		partition(indices, pivot + 1, i1);

	}


	int test() {

		int values[12];
		for (int k = 0; k < 20; k++) {
			for (int i = 0; i < 12; i++) {
				values[i] = rand() % 100;
			}
			partition(values, 0, 12);
		}
		return 0;
	}*/

	void paintEvent(wxPaintEvent& evt)
	{
	   //test();
		wxPaintDC dc(this);
		render(dc);
	}

	inline double fastPow(double a, double b) {
		union {
			double d;
			int x[2];
		} u = { a };
		u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
		u.x[0] = 0;
		return u.d;
	}

	void delete_object(wxCommandEvent& event) {
		if (selected_object < 0) return;
		if (selected_object >= raytracer.s.objects.size()) return;

		stop_render();
		raytracer.s.deleteObject(selected_object);
		start_render();
	}

	void colorAnisotropy(wxCommandEvent& event) {
		if (selected_object < 0) return;
		if (selected_object >= raytracer.s.objects.size()) return;

		stop_render();
		raytracer.s.objects[selected_object]->colorAnisotropy();
		start_render();
	}

	void randomColors(wxCommandEvent& event) {
		if (selected_object < 0) return;
		if (selected_object >= raytracer.s.objects.size()) return;

		stop_render();
		raytracer.s.objects[selected_object]->randomColors();
		start_render();
	}

	void add_keyframe(wxCommandEvent& event) {
		if (selected_object < 0) return;
		if (selected_object >= raytracer.s.objects.size()) return;

		raytracer.s.objects[selected_object]->add_keyframe(raytracer.s.current_frame);
	}

	void launch_render(wxCommandEvent& event) {
		stop_render();
		PerfChrono perf;
		perf.Start();
		raytracer.stopped = false;
		raytracer.render_image();
		raytracer.stopped = true;
		double t = perf.GetDiffMs();
		FILE* f = fopen("time.txt", "a+");
		fprintf(f, "%e\n", t);
		fclose(f);
		start_render();
	}

	void render_video(wxCommandEvent& event) {
		stop_render();
		PerfChrono perf;
		perf.Start();
		raytracer.stopped = false;
		for (int i = 0; i < raytracer.s.nbframes; i++) {
			raytracer.clear_image();
			raytracer.s.current_frame = i;
			raytracer.s.current_time = i*raytracer.s.duration / (double)raytracer.s.nbframes;
			raytracer.render_image();
		}
		raytracer.stopped = true;
		start_render();
	}

	void update_textures_and_render(wxCommandEvent& event);
	void update_parameters_and_render(wxCommandEvent& event);
	void update_gui();

	std::vector<double> extrapolated_image;
	std::vector<bool> computed, computed2;

	void render(wxDC& dc);



	void paintNow()
	{
		wxClientDC dc(this);
		render(dc);
	}
	void stop_render() {
		raytracer.stopRender();		
		//WaitForSingleObject(compute_thread, INFINITE);
		compute_thread2.join();
	}

	void start_render() {
		//compute_thread = CreateThread(NULL, 0, doRender, (void*)this, 0, NULL);
		//SetThreadPriority(compute_thread, THREAD_PRIORITY_HIGHEST);
		raytracer.clear_image();
		raytracer.stopped = false;
		compute_thread2 = std::thread(doRender, (void*)this);		
	}

	void mouseDown(wxMouseEvent& event) {
		mouse_init_x = event.m_x; // coordinates before dragging the mouse
		mouse_init_y = event.m_y;
	}

	void mouseMoved(wxMouseEvent& event)
	{
		mouse_x = event.m_x;
		mouse_y = event.m_y;

		if (mouse_x > displayW + 10) return;
		if (mouse_y > displayH + 10) return;

		if ((event.LeftIsDown() || event.RightIsDown() || event.MiddleIsDown()) && event.Dragging())
		{
			left_mouse_down = event.LeftIsDown();
			right_mouse_down = event.RightIsDown();
			middle_mouse_down = event.MiddleIsDown();
			alt_down = event.AltDown();
			shift_down = event.ShiftDown();
			wasDragging = true;

			int dx = mouse_x - mouse_prev_x;
			int dy = -mouse_y + mouse_prev_y;
			
			if (shift_down) {
				if (left_mouse_down) {
					if (selected_object >= 0 && selected_object < raytracer.s.objects.size()) {
						//raytracer.s.objects[selected_object]->max_translation += Vector(dx, dy, 0.)*0.3;
						Vector camera_right = cross(raytracer.cam.direction, raytracer.cam.up);
						raytracer.s.objects[selected_object]->max_translation += dx*selected_object_t*tan(raytracer.cam.fov/ displayW)*camera_right + dy*selected_object_t*tan(raytracer.cam.fov / displayH)*raytracer.cam.up;
					}
				}
				if (middle_mouse_down) {
					if (selected_object >= 0 && selected_object < raytracer.s.objects.size()) {
						//raytracer.s.objects[selected_object]->max_rotation += Vector(0., 0., dy*M_PI / 180.*3.);
						Matrix<3, 3> R = createRotationMatrixZ(dy*M_PI / 180.*3.);

						raytracer.s.objects[selected_object]->mat_rotation = R * raytracer.s.objects[selected_object]->mat_rotation;
					}
				}
				if (right_mouse_down) {
					if (selected_object >= 0 && selected_object < raytracer.s.objects.size()) {
						//raytracer.s.objects[selected_object]->max_rotation += Vector(dy*M_PI / 180.*3., dx*M_PI / 180.*3., 0.);
						Matrix<3, 3> R = createRotationMatrixX(dy*M_PI / 180.*3.)*createRotationMatrixY(dx*M_PI / 180.*3.);
						raytracer.s.objects[selected_object]->mat_rotation = R * raytracer.s.objects[selected_object]->mat_rotation;
					}
				}
			} else { // camera motion
				if (left_mouse_down) {

					if (abs(mouse_x - displayW) < 10 && abs(mouse_y - displayH) < 10) {
						displayW = mouse_x;
						displayH = mouse_y;
					} else
						raytracer.cam.rotateAroundAxes(dx*M_PI / 180 * 1, dy*M_PI / 180 * 1, 1);
				}
				if (right_mouse_down) {
					Vector camera_right = cross(raytracer.cam.direction, raytracer.cam.up);
					raytracer.cam.position += dx*camera_right + dy*raytracer.cam.up;

				}
			}			
			stop_render();
			start_render();
										
		} else {
			if (!event.Dragging()) {

			}
		}
		mouse_prev_x = mouse_x;
		mouse_prev_y = mouse_y;

		paintNow();
	}

	void mouseUp(wxMouseEvent& event) {
				
		if ( (event.m_x > displayW + 10) || (event.m_y > displayW + 10)) {
			wasDragging = false;
			return;
		}	
		

		/*if (!wasDragging)*/ {
			Ray r = raytracer.cam.generateDirection((displayH - (mouse_init_y - 1.))*(double)raytracer.H / displayH, mouse_init_x*(double)raytracer.W / displayW, 0, 0, 0, 0, 0, raytracer.W, raytracer.H);

			Vector P;
			int new_selected, new_tri;
			double new_t;
			MaterialValues newmat;
			bool has_inter = raytracer.s.intersection(r, P, new_selected, new_t, newmat, new_tri);
			if (has_inter) {
				selected_object = new_selected;
				selected_object_t = new_t;
				selected_tri = new_tri;
				update_gui();
			}			
			mouse_prev_x = mouse_x;
			mouse_prev_y = mouse_y;

		}
		
		left_mouse_down = false;
		if (!event.AltDown())
			alt_down = false;

		wasDragging = false;
	}

	void mouseWheel(wxMouseEvent& event) {
		if (event.ShiftDown()) {
			if (selected_object >= 0 && selected_object < raytracer.s.objects.size()) {
				raytracer.s.objects[selected_object]->scale *= (event.GetWheelRotation() > 0) ? 1.1 : (1 / 1.1);
				stop_render();
				start_render();
			}
		} else {
			double dist =  (event.GetWheelRotation() > 0 ? 1 : -1)*2.;
			raytracer.cam.position += dist*raytracer.cam.direction;
			stop_render();
			start_render();
		}
		paintNow();
	}

	bool should_terminate;
	bool left_mouse_down, right_mouse_down, middle_mouse_down, alt_down, shift_down, wasDragging;
	int mouse_x, mouse_y, mouse_init_x, mouse_init_y, mouse_prev_x, mouse_prev_y, selected_object, selected_tri;
	double selected_object_t;
	std::vector<unsigned char> cur_img;
	int W, H;
	int displayW, displayH;
	//HANDLE compute_thread;
	std::thread compute_thread2;

	RaytracerApp* raytracer_app;
	Raytracer raytracer;
	wxImage screenImage;
	wxBitmap bmpBuf;
	std::vector<unsigned char> gamma_corrected_image;

	DECLARE_EVENT_TABLE()
};

static void doRender(void* Param) {
	Raytracer* raytracer = &(((RenderPanel*)Param)->raytracer);
	raytracer->render_image();
}

BEGIN_EVENT_TABLE(RenderPanel, wxPanel)
EVT_PAINT(RenderPanel::paintEvent)
EVT_MOTION(RenderPanel::mouseMoved)
EVT_LEFT_DOWN(RenderPanel::mouseDown)
EVT_RIGHT_DOWN(RenderPanel::mouseDown)
EVT_LEFT_UP(RenderPanel::mouseUp)
EVT_RIGHT_UP(RenderPanel::mouseUp)
EVT_MOUSEWHEEL(RenderPanel::mouseWheel)
END_EVENT_TABLE()


class RaytracerApp : public wxApp
{
public:
	virtual bool OnInit() ;
	RenderPanel* renderPanel;
	wxGauge *progressBar;
	wxBookCtrlBase *m_bookCtrl;
	wxCheckBox *show_edges, *interp_normals, *transparent, *flipnormals, *isLenticularCheck, *ghost;
	wxTextCtrl *objectName, *envmapName, *backgroundName;
	wxSlider *fov_slider, *aperture_slider, /**ks_slider,*/ *filter_slider, *fogdensity_slider, *envmapintensity_slider, *lightintensity_slider, *focus_slider, *maxangle_slider, *time_slider;
	wxListCtrl *m_AlbedoFile, *m_SpecularFile, *m_NormalFile, *m_AlphaFile, *m_RoughnessFile, *m_TranspFile, *m_RefrFile;
	//wxColourPickerCtrl *albedoColorPicker;
	wxRadioButton *uniformFogRadio, *expFogRadio;
	wxColourDialog* colPicker;
	wxFileDialog* texOpenDlg;
	wxSpinCtrl *bounces, *renderwidth, *renderheight, *nbrays, *nbviews, *nbpixslice, *nbframesctrl;
	wxSpinCtrlDouble *duration;
	wxSpinCtrlDouble* refractionIndex;
	wxButton *deleteObject, *launchRender, *addKeyframe, *renderVideo, *colorAnisotropy, *randomColors;
	wxToggleButton*recordKeyframes;
	wxStaticText* infoModel;
	wxStaticText* infoPerf;
	bool render_loop_on;

	void activateRenderLoop(bool on);

	void OnEraseBackGround(wxEraseEvent& event) {
		//wxDC * TheDC = event.GetDC();		
		event.Skip(); // don't call -- not needed, we already erased the background
	}

	void OnIdle(wxIdleEvent& evt)
	{
		if (render_loop_on)
		{
			renderPanel->paintNow();
			evt.RequestMore(); // render continuously, not only once on idle				
		}
	}
	void OnUnhandledException() {
		throw;
	};

	DECLARE_EVENT_TABLE()
};

BEGIN_EVENT_TABLE(RaytracerApp, wxApp)
EVT_IDLE(RaytracerApp::OnIdle)
END_EVENT_TABLE()



