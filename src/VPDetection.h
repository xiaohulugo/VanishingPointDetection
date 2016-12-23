#ifndef _VP_DETECTION_H_
#define _VP_DETECTION_H_
#pragma once

#include <stdio.h>
#include <opencv2/opencv.hpp> 

struct LineInfo
{
	cv::Mat_<double> para;
	double length;
	double orientation;
};

class VPDetection
{
public:
	VPDetection(void);
	~VPDetection(void);

	void run( std::vector<std::vector<double> > &lines, cv::Point2d pp, double f, std::vector<cv::Point3d> &vps, std::vector<std::vector<int> > &clusters );

	void getVPHypVia2Lines( std::vector<std::vector<cv::Point3d> >  &vpHypo );

	void getSphereGrids( std::vector<std::vector<double> > &sphereGrid );

	void getBestVpsHyp( std::vector<std::vector<double> > &sphereGrid, std::vector<std::vector<cv::Point3d> >  &vpHypo, std::vector<cv::Point3d> &vps  );

	void lines2Vps( double thAngle, std::vector<cv::Point3d> &vps, std::vector<std::vector<int> > &clusters );

private:
	std::vector<std::vector<double> > lines;
	std::vector<LineInfo> lineInfos;
	cv::Point2d pp;
	double f;
	double noiseRatio;
};

#endif // _VP_DETECTION_H_
