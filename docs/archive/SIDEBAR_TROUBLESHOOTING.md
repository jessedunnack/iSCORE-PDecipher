# Sidebar Toggle Troubleshooting Guide

## Current Implementation Status

### ✅ Completed Fixes
1. **Sticky Note Content Enhanced**:
   - Now clearly describes evaluation of gene mutation/perturbation effects
   - Mentions analysis at both DEG and functional enrichment levels
   - Explicitly states that sidebar is collapsible via hamburger menu (☰)

2. **Sidebar Click Handler Improved**:
   - Added `stopPropagation()` to prevent event conflicts
   - Disabled auto-management on Overview page
   - Enhanced debugging with console logs

## If Clicking Still Doesn't Work

### Quick Debugging Steps

1. **Open Browser Console** (F12) and check for:
   - Look for `[SIDEBAR]` messages when clicking
   - Any JavaScript errors in red

2. **Try Force Toggle via Console**:
   ```javascript
   // Paste this in console to manually toggle
   toggleSidebar();
   ```

3. **Check Button Visibility**:
   ```javascript
   // Check if button exists and is visible
   $('.sidebar-toggle').length  // Should return 1
   $('.sidebar-toggle').is(':visible')  // Should return true
   ```

### Potential Issues & Solutions

#### Z-Index Conflict
If the sticky note is overlapping the toggle button:
```css
/* Add to custom CSS */
.sidebar-toggle {
  z-index: 1051 !important;  /* Higher than sticky note */
}
```

#### Event Handler Not Attached
If click handler isn't working:
```javascript
// Force re-attach handler
$('.sidebar-toggle').off('click').on('click', function(e) {
  e.preventDefault();
  toggleSidebar();
});
```

#### CSS Pointer Events
If button appears but can't be clicked:
```css
.sidebar-toggle {
  pointer-events: auto !important;
}
```

## Alternative Workarounds

### 1. Keyboard Shortcut
Add this to enable Ctrl+B to toggle sidebar:
```javascript
$(document).keydown(function(e) {
  if (e.ctrlKey && e.which === 66) {  // Ctrl+B
    e.preventDefault();
    toggleSidebar();
  }
});
```

### 2. Direct Function Call
Create a temporary button in console:
```javascript
$('body').append('<button onclick="toggleSidebar()" style="position:fixed;top:10px;right:10px;z-index:9999">Toggle</button>');
```

## Expected Behavior

1. **On Overview Page**:
   - Sticky note appears on left
   - Sidebar starts collapsed
   - Manual clicking should always work
   - No auto-expansion/collapse

2. **On Other Pages**:
   - Sidebar auto-expands when needed
   - Manual override still works
   - State persists for 1 second after manual toggle

## Console Commands for Testing

```javascript
// Check current sidebar state
console.log('Sidebar collapsed:', sidebarState.isCollapsed);

// Check if on Overview page
console.log('Current tab:', getCurrentTab());

// Force expand
toggleSidebar(true);

// Force collapse  
toggleSidebar(false);

// Check if button click is registered
$('.sidebar-toggle').click(function() { 
  console.log('Button clicked!'); 
});
```

## If All Else Fails

The sidebar functionality is not critical for using the app. The global settings affect all visualizations, so you can:
1. Set them once when expanded
2. Navigate to different analysis tabs
3. The settings will persist across all tabs

The sticky note provides the key information about using these settings, which is the most important user guidance feature.