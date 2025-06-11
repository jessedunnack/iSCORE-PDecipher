// JavaScript enhancements for PD Enrichment Explorer

// Add smooth scrolling
$(document).on('click', 'a[href^="#"]', function (event) {
    event.preventDefault();
    $('html, body').animate({
        scrollTop: $($.attr(this, 'href')).offset().top
    }, 500);
});

// Add tooltips
$(document).ready(function(){
    $('[data-toggle="tooltip"]').tooltip();
});

// Custom notification function
function showNotification(message, type = 'info') {
    Shiny.notifications.show({
        html: message,
        type: type,
        duration: 3000
    });
}

// Add keyboard shortcuts
$(document).keydown(function(e) {
    // Ctrl+S to save/export (when in export tab)
    if (e.ctrlKey && e.which === 83) {
        e.preventDefault();
        if ($('#sidebar_menu .active').attr('data-value') === 'export') {
            $('#export-export_data').click();
        }
    }
});

// Add loading overlay functions
function showLoading() {
    $('body').append('<div class="loading-overlay"><div class="spinner-border text-primary" role="status"><span class="sr-only">Loading...</span></div></div>');
}

function hideLoading() {
    $('.loading-overlay').remove();
}

// Enhance file upload area with drag and drop
$(document).on('dragover', '.file-upload-area', function(e) {
    e.preventDefault();
    $(this).addClass('drag-over');
});

$(document).on('dragleave', '.file-upload-area', function(e) {
    e.preventDefault();
    $(this).removeClass('drag-over');
});

// Add console message
console.log('%cPD Enrichment Explorer', 'color: #3c8dbc; font-size: 20px; font-weight: bold;');
console.log('Version 1.0.0 - iSCORE-PDecipher Project');