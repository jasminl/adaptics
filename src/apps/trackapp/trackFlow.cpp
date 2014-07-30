#include "trackFlow.h"

using namespace std;

trackFlow::trackFlow(trackFilter* filt, trackMatch* match):m_currentX(0),m_currentY(0),m_currentScale(0)
{
	m_featFilter = filt;
	m_targetFilter = filt;	//TODO: why is this the same as previous?
	m_match = match;
}
