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

// ===== SIDEBAR TOGGLE FUNCTIONALITY =====

// Sidebar state management
let sidebarState = {
    isCollapsed: true,  // Start collapsed by default 
    userManualToggle: false,  // Track if user manually toggled
    lastAutoState: null  // Track last automatic state to avoid conflicts
};

// Function to get current active tab
function getCurrentTab() {
    // Try multiple methods to detect the active tab
    let mainTab = null;
    
    // Method 1: Check for active nav-link with data-value
    const activeNavLink = $('#main_sections .nav-link.active');
    if (activeNavLink.length > 0) {
        mainTab = activeNavLink.attr('data-value') || activeNavLink.attr('href')?.replace('#', '');
    }
    
    // Method 2: Check for active tab-pane
    if (!mainTab) {
        const activePane = $('#main_sections .tab-pane.active');
        if (activePane.length > 0) {
            mainTab = activePane.attr('data-value') || activePane.attr('id');
        }
    }
    
    // Method 3: Check Shiny's internal state if available
    if (!mainTab && window.Shiny && window.Shiny.inputBindings) {
        const tabInput = window.Shiny.shinyapp && window.Shiny.shinyapp.$inputValues && 
                        window.Shiny.shinyapp.$inputValues['main_sections'];
        if (tabInput) {
            mainTab = tabInput;
        }
    }
    
    // Method 4: Parse from URL hash
    if (!mainTab) {
        const hash = window.location.hash;
        if (hash && hash.length > 1) {
            mainTab = hash.substring(1);
        }
    }
    
    // Default fallback
    mainTab = mainTab || 'overview';
    
    // Clean up the tab value
    mainTab = mainTab.replace(/^#/, '').replace(/^tab-/, '');
    
    return mainTab;
}

// Function to determine if current page needs settings expanded
function shouldExpandSettings(tabValue) {
    const settingsRequiredTabs = [
        'de_genes_section',      // All DE Genes subtabs need settings
        'enrichment_section',    // All Functional Enrichment subtabs need settings  
        'signature_nomination'   // All Signature Nomination subtabs need settings
    ];
    
    return settingsRequiredTabs.includes(tabValue);
}

// Function to toggle sidebar state
function toggleSidebar(forceState = null, isAutomatic = false) {
    const sidebar = $('.sidebar-fixed');
    const mainContent = $('.main-content');
    
    let newState;
    if (forceState !== null) {
        newState = forceState;  // true = expanded, false = collapsed
    } else {
        newState = !sidebarState.isCollapsed;  // Toggle current state
    }
    
    // Update visual state
    if (newState) {
        // Expand sidebar
        sidebar.removeClass('collapsed');
        mainContent.removeClass('expanded');
        sidebarState.isCollapsed = false;
        console.log('[SIDEBAR] Expanded' + (isAutomatic ? ' (automatic)' : ' (manual)'));
    } else {
        // Collapse sidebar  
        sidebar.addClass('collapsed');
        mainContent.addClass('expanded');
        sidebarState.isCollapsed = true;
        console.log('[SIDEBAR] Collapsed' + (isAutomatic ? ' (automatic)' : ' (manual)'));
    }
    
    // Track if this was a manual user action
    if (!isAutomatic) {
        sidebarState.userManualToggle = true;
        setTimeout(() => { sidebarState.userManualToggle = false; }, 1000);  // Reset after 1 second
    }
    
    sidebarState.lastAutoState = isAutomatic ? newState : sidebarState.lastAutoState;
}

// Function to manage sidebar based on current tab
function manageSidebarForTab() {
    const currentTab = getCurrentTab();
    const shouldExpand = shouldExpandSettings(currentTab);
    
    console.log('[SIDEBAR] Current tab:', currentTab, '| Should expand:', shouldExpand);
    
    // Only auto-manage if user hasn't manually toggled recently
    if (!sidebarState.userManualToggle) {
        const needsChange = (shouldExpand && sidebarState.isCollapsed) || 
                           (!shouldExpand && !sidebarState.isCollapsed);
                           
        if (needsChange) {
            toggleSidebar(shouldExpand, true);  // true = automatic
        }
    }
}

// Initialize sidebar toggle functionality
$(document).ready(function() {
    // Set initial collapsed state
    toggleSidebar(false, true);  // Start collapsed
    
    // Manual toggle button click handler
    $(document).on('click', '.sidebar-toggle', function(e) {
        e.preventDefault();
        console.log('[SIDEBAR] Manual toggle clicked');
        toggleSidebar();  // Manual toggle
    });
    
    // Monitor main tab changes
    $(document).on('shown.bs.tab', 'a[data-toggle="tab"]', function (e) {
        setTimeout(manageSidebarForTab, 100);  // Small delay to ensure tab is fully loaded
    });
    
    // Monitor radio button changes for main sections
    $(document).on('change', 'input[name="main_sections"]', function() {
        setTimeout(manageSidebarForTab, 100);
    });
    
    // Also monitor clicks on tab links directly
    $(document).on('click', '.nav-tabs .nav-link, .nav-pills .nav-link', function() {
        setTimeout(manageSidebarForTab, 200);
    });
    
    // Monitor for dynamic tab updates via Shiny
    if (window.Shiny) {
        // Hook into Shiny's input change system
        $(document).on('shiny:inputchanged', function(event) {
            if (event.name === 'main_sections') {
                console.log('[SIDEBAR] Shiny tab changed to:', event.value);
                setTimeout(manageSidebarForTab, 150);
            }
        });
        
        // Custom message handler for manual tab updates
        Shiny.addCustomMessageHandler('updateMainTab', function(message) {
            setTimeout(manageSidebarForTab, 150);
        });
        
        // Monitor when Shiny is ready
        $(document).on('shiny:sessioninitialized', function() {
            console.log('[SIDEBAR] Shiny session initialized');
            setTimeout(manageSidebarForTab, 500);  // Check initial state
        });
    }
    
    // Enhanced periodic check for tab changes (more responsive)
    let lastCheckedTab = null;
    setInterval(function() {
        const currentTab = getCurrentTab();
        if (currentTab !== lastCheckedTab) {
            lastCheckedTab = currentTab;
            console.log('[SIDEBAR] Tab change detected via polling:', currentTab);
            manageSidebarForTab();
        }
    }, 1000);  // Check every second
    
    console.log('[SIDEBAR] Auto-management initialized');
});

// Add console message
console.log('%cPD Enrichment Explorer', 'color: #3c8dbc; font-size: 20px; font-weight: bold;');
console.log('Version 1.0.0 - iSCORE-PDecipher Project');
console.log('%cSidebar auto-management enabled', 'color: #28a745; font-weight: bold;');