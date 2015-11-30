jQuery(document).ready(function() {

	viewer = new Viewer('#layer_list', '.layer_settings');
	viewer.addView('#view_axial', Viewer.AXIAL);
	viewer.addView('#view_coronal', Viewer.CORONAL);
	viewer.addView('#view_sagittal', Viewer.SAGITTAL);
	//viewer.addSlider('opacity', '.slider#opacity', 'horizontal', 0, 1, 1, 0.05);
	//viewer.addSlider('pos-threshold', '.slider#pos-threshold', 'horizontal', 0, 1, 0, 0.01);
	//viewer.addSlider('neg-threshold', '.slider#neg-threshold', 'horizontal', 0, 1, 0, 0.01);
	//viewer.addSlider("nav-xaxis", ".slider#nav-xaxis", "horizontal", 0, 1, 0.5, 0.01, Viewer.XAXIS);
	//viewer.addSlider("nav-yaxis", ".slider#nav-yaxis", "vertical", 0, 1, 0.5, 0.01, Viewer.YAXIS);
	//viewer.addSlider("nav-zaxis", ".slider#nav-zaxis", "vertical", 0, 1, 0.5, 0.01, Viewer.ZAXIS);

	viewer.addDataField('voxelValue', '#data_current_value')
	viewer.clear()   // Paint canvas background while images load
	images = [
		{
			'url': 'data/flip_crop.nii.gz',
			'name': 'Mean n592',
			'colorPalette': 'grayscale',
			'cache': false,
			'intent': 'Intensity:'
		},
		{
			'url': 'data/AT_mean_T1ToD_NORSRS_0025_cluster.nii.gz',
			'name': 'Neg',
			'sign': 'both',
			'colorPalette': 'red-blue',
			'cache': false,
			'intent': 'Intensity:'
		}
	]
	viewer.loadImages(images);

});
