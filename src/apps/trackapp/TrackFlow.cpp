#include "TrackFlow.h"

using namespace std;

TrackFlow::TrackFlow(trackFilter* filt, TrackMatch* match)
: _cur_x(0), _cur_y(0), _cur_scale(0)
{
	_feat_filter = filt;
	_target_filter = filt;	//TODO: why is this the same as previous?
	_match = match;
}
