# ----------------------------------
#  Upload and share files with google drive
# ----------------------------------

#' Upload and share
#'
#' Upload and share (view only) a file on google drive.  Will authenticate in
#' your browser the first time you use it
#'
#' @param file local file to upload
#' @param path path on google drive to upload to
#' @param remove_existing logical, remove any existing files with the same name in the given path
#'
#' @return a list with the viewing and download link
#'
#' @export
upload_and_share = function(file, path, remove_existing = TRUE) {
	dir = googledrive::drive_get(path)
	if (nrow(dir) == 0) {
		message("upload dir doesn't exist, creating")
		dir = googledrive::drive_mkdir(path)
	}
	existing = googledrive::drive_ls(dir) %>%
		dplyr::filter(name == file)
	if (nrow(existing) != 0 & remove_existing) {
		message("Overwriting existing file...")
		googledrive::drive_trash(existing)
	} else if (nrow(existing) != 0 & !remove_existing) {
		stop("File existing and remove_existing is FALSE")
	}
	upload = googledrive::drive_upload(file, dir)
	shared = googledrive::drive_share(upload, role = "reader", type = "anyone")
	return(list(
		viewLink = purrr::pluck(shared$drive_resource, 1, "webViewLink"),
		downloadLink = purrr::pluck(shared$drive_resource, 1, "webContentLink")

	))
}
